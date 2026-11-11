"""Offline tests for the SQLite layer. Run: pytest -q"""
import os

import pytest

import store_sqlite as db
from build_dashboard import load_records

SAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_data", "pubmed_sample.xml")


def rec(pmid, title="A paper", doi="", year=2023, abstract="", authors=("A B", "C D"), mesh=()):
    return {"pmid": pmid, "title": title, "abstract": abstract, "journal": "J", "year": year,
            "doi": doi, "authors": list(authors), "mesh": list(mesh), "pub_types": [], "url": ""}


@pytest.fixture
def conn():
    c = db.connect(":memory:")
    yield c
    c.close()


def test_load_sample_dedupes_and_counts(conn):
    r = db.add_records(conn, load_records(SAMPLE), "sample")
    assert r["fetched"] == 60 and r["new"] == 60 and r["duplicates_in_batch"] == 3
    assert db.stats(conn)["papers"] == 60


def test_rerun_adds_nothing_new(conn):
    recs = load_records(SAMPLE)
    db.add_records(conn, recs, "q")
    again = db.add_records(conn, recs, "q")
    assert again["new"] == 0 and again["already_stored"] == 60
    assert db.stats(conn)["papers"] == 60 and db.stats(conn)["searches"] == 2


def test_same_doi_different_pmid_is_not_stored_twice(conn):
    db.add_records(conn, [rec("1", doi="10.1/X")])
    r = db.add_records(conn, [rec("2", doi="10.1/x")])         # DOI differs only by case
    assert r["new"] == 0 and db.stats(conn)["papers"] == 1


def test_missing_doi_does_not_collide(conn):
    db.add_records(conn, [rec("1"), rec("2")])                 # both without a DOI
    assert db.stats(conn)["papers"] == 2


def test_authors_shared_between_papers_stored_once(conn):
    db.add_records(conn, [rec("1", authors=("Kim J", "Lee S")), rec("2", authors=("Kim J",))])
    assert db.stats(conn)["authors"] == 2
    order = [r[0] for r in conn.execute(
        "SELECT a.name FROM paper_authors pa JOIN authors a ON a.id = pa.author_id "
        "JOIN papers p ON p.id = pa.paper_id WHERE p.pmid = '1' ORDER BY pa.position")]
    assert order == ["Kim J", "Lee S"]                         # author order preserved


def test_topic_tags_and_search(conn):
    topics = {"Safety": ["hepatotox"], "ML": ["random forest"]}
    db.add_records(conn, [rec("1", "Hepatotoxicity signals", year=2021),
                          rec("2", "Random forest model", year=2024),
                          rec("3", "Unrelated study", year=2022)], topics=topics)
    assert db.stats(conn)["by_topic"] == {"ML": 1, "Safety": 1}
    assert [r["pmid"] for r in db.search(conn, topic="Safety")] == ["1"]
    assert [r["pmid"] for r in db.search(conn, text="MODEL")] == ["2"]      # case-insensitive
    assert [r["pmid"] for r in db.search(conn, year_from=2022)] == ["2", "3"]   # newest first


def test_search_is_safe_against_sql_injection(conn):
    db.add_records(conn, [rec("1", "Normal title")])
    assert db.search(conn, text="'; DROP TABLE papers; --") == []
    assert db.stats(conn)["papers"] == 1                       # table still there


def test_cli_roundtrip(tmp_path, capsys):
    path = str(tmp_path / "t.db")
    assert db.main(["load", SAMPLE, "--db", path]) == 0
    assert db.main(["stats", "--db", path]) == 0
    out = capsys.readouterr().out
    assert "'new': 60" in out and "'papers': 60" in out
