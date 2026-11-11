"""
SQLite storage layer for the literature pipeline.

    python store_sqlite.py load  sample_data/pubmed_sample.xml --db literature.db --query "demo search"
    python store_sqlite.py search --db literature.db --text hepatotox --topic "Machine Learning"
    python store_sqlite.py stats --db literature.db

Why a database: records survive between runs, re-running a search only adds what is new
(PMID and DOI are unique), and the dashboard or any other tool can query one file.

Schema
    papers(id, pmid UNIQUE, doi UNIQUE-when-present, title, abstract, journal, year, url, added_at)
    authors(id, name UNIQUE)              paper_authors(paper_id, author_id, position)
    tags(id, name UNIQUE)                 paper_tags(paper_id, tag_id)
    searches(id, query, run_at, n_fetched, n_new)

Standard library plus the project's own modules (no server needed).
"""
from __future__ import annotations

import argparse
import sqlite3
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional

from analyze_literature import DEFAULT_TOPICS, _compile, dedupe, tag_topics
from build_dashboard import load_records

SCHEMA = """
CREATE TABLE IF NOT EXISTS papers (
    id        INTEGER PRIMARY KEY,
    pmid      TEXT UNIQUE,
    doi       TEXT,
    title     TEXT NOT NULL,
    abstract  TEXT NOT NULL DEFAULT '',
    journal   TEXT NOT NULL DEFAULT '',
    year      INTEGER,
    url       TEXT NOT NULL DEFAULT '',
    added_at  TEXT NOT NULL
);
CREATE UNIQUE INDEX IF NOT EXISTS ux_papers_doi ON papers(doi) WHERE doi IS NOT NULL;
CREATE TABLE IF NOT EXISTS authors (id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL);
CREATE TABLE IF NOT EXISTS paper_authors (
    paper_id  INTEGER NOT NULL REFERENCES papers(id) ON DELETE CASCADE,
    author_id INTEGER NOT NULL REFERENCES authors(id),
    position  INTEGER NOT NULL,
    PRIMARY KEY (paper_id, author_id)
);
CREATE TABLE IF NOT EXISTS tags (id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL);
CREATE TABLE IF NOT EXISTS paper_tags (
    paper_id INTEGER NOT NULL REFERENCES papers(id) ON DELETE CASCADE,
    tag_id   INTEGER NOT NULL REFERENCES tags(id),
    PRIMARY KEY (paper_id, tag_id)
);
CREATE TABLE IF NOT EXISTS searches (
    id INTEGER PRIMARY KEY, query TEXT NOT NULL, run_at TEXT NOT NULL,
    n_fetched INTEGER NOT NULL, n_new INTEGER NOT NULL
);
"""


def connect(path: str = ":memory:") -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA)
    return conn


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _exists(conn, rec: dict) -> bool:
    pmid, doi = rec.get("pmid") or None, (rec.get("doi") or "").strip().lower() or None
    row = conn.execute("SELECT 1 FROM papers WHERE (pmid IS NOT NULL AND pmid = ?) "
                       "OR (doi IS NOT NULL AND doi = ?)", (pmid, doi)).fetchone()
    return row is not None


def _get_id(conn, table: str, name: str) -> int:
    conn.execute(f"INSERT OR IGNORE INTO {table}(name) VALUES (?)", (name,))
    return conn.execute(f"SELECT id FROM {table} WHERE name = ?", (name,)).fetchone()[0]


def add_records(conn: sqlite3.Connection, records: Iterable[dict], query: str = "",
                topics: Optional[Dict[str, List[str]]] = None) -> Dict[str, int]:
    """Insert records that are not already stored (by PMID or DOI). Returns counts."""
    patterns = _compile(topics or DEFAULT_TOPICS)
    unique, n_dup = dedupe(list(records))
    n_new = 0
    with conn:                                              # one transaction
        for rec in unique:
            if _exists(conn, rec):
                continue
            cur = conn.execute(
                "INSERT INTO papers(pmid, doi, title, abstract, journal, year, url, added_at) "
                "VALUES (?,?,?,?,?,?,?,?)",
                (rec.get("pmid") or None, (rec.get("doi") or "").strip().lower() or None,
                 rec.get("title", ""), rec.get("abstract", ""), rec.get("journal", ""),
                 rec.get("year") or None, rec.get("url", ""), _now()))
            pid = cur.lastrowid
            for pos, name in enumerate(rec.get("authors", [])):
                conn.execute("INSERT OR IGNORE INTO paper_authors VALUES (?,?,?)",
                             (pid, _get_id(conn, "authors", name), pos))
            for tag in tag_topics(rec, patterns):
                conn.execute("INSERT OR IGNORE INTO paper_tags VALUES (?,?)",
                             (pid, _get_id(conn, "tags", tag)))
            n_new += 1
        n_in = len(unique)
        conn.execute("INSERT INTO searches(query, run_at, n_fetched, n_new) VALUES (?,?,?,?)",
                     (query, _now(), n_in, n_new))
    return {"fetched": n_in, "new": n_new, "already_stored": n_in - n_new, "duplicates_in_batch": n_dup}


def search(conn: sqlite3.Connection, text: str = "", topic: str = "", year_from: Optional[int] = None,
           limit: int = 50) -> List[sqlite3.Row]:
    """Papers whose title/abstract contain `text`, optionally in a topic and from a year on.
    All values are passed as bound parameters, never pasted into the SQL."""
    sql = ["SELECT DISTINCT p.* FROM papers p"]
    args: list = []
    if topic:
        sql.append("JOIN paper_tags pt ON pt.paper_id = p.id JOIN tags t ON t.id = pt.tag_id AND t.name = ?")
        args.append(topic)
    sql.append("WHERE 1=1")
    if text:
        like = f"%{text.lower()}%"
        sql.append("AND (lower(p.title) LIKE ? OR lower(p.abstract) LIKE ?)")
        args += [like, like]
    if year_from:
        sql.append("AND p.year >= ?")
        args.append(year_from)
    sql.append("ORDER BY p.year DESC, p.id LIMIT ?")
    args.append(limit)
    return conn.execute(" ".join(sql), args).fetchall()


def stats(conn: sqlite3.Connection) -> dict:
    one = lambda q: conn.execute(q).fetchone()[0]           # noqa: E731
    return {
        "papers": one("SELECT COUNT(*) FROM papers"),
        "authors": one("SELECT COUNT(*) FROM authors"),
        "by_year": {r["year"]: r["n"] for r in conn.execute(
            "SELECT year, COUNT(*) n FROM papers WHERE year IS NOT NULL GROUP BY year ORDER BY year")},
        "by_topic": {r["name"]: r["n"] for r in conn.execute(
            "SELECT t.name, COUNT(*) n FROM tags t JOIN paper_tags pt ON pt.tag_id = t.id "
            "GROUP BY t.name ORDER BY n DESC, t.name")},
        "searches": one("SELECT COUNT(*) FROM searches"),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="SQLite store for literature records.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("load", help="load an XML/JSON records file")
    p.add_argument("input"); p.add_argument("--db", default="literature.db"); p.add_argument("--query", default="")
    p = sub.add_parser("search", help="search stored papers")
    p.add_argument("--db", default="literature.db"); p.add_argument("--text", default="")
    p.add_argument("--topic", default=""); p.add_argument("--year-from", type=int); p.add_argument("--limit", type=int, default=20)
    p = sub.add_parser("stats", help="counts by year and topic")
    p.add_argument("--db", default="literature.db")
    args = ap.parse_args(argv)

    conn = connect(args.db)
    if args.cmd == "load":
        print(add_records(conn, load_records(args.input), args.query or args.input))
    elif args.cmd == "search":
        for r in search(conn, args.text, args.topic, args.year_from, args.limit):
            print(f"{r['year']}  {r['pmid']}  {r['title']}")
    else:
        print(stats(conn))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
