"""
Clean, de-duplicate, tag and summarise literature records.

    from analyze_literature import dedupe, to_dataframe, summarize

Steps
  dedupe()        remove repeats (same DOI, else same PMID, else same normalised title)
  to_dataframe()  one row per paper + topic tags from a keyword dictionary
  summarize()     counts per year / journal / topic / MeSH term for the dashboard

Topic tagging is plain keyword matching on title + abstract + MeSH terms. It is transparent
and easy to audit, but it is a first-pass screen: a person should review what it flags.
"""
from __future__ import annotations

import json
import re
from typing import Dict, List, Tuple

import pandas as pd

# keyword prefixes: "hepatotox" also matches "hepatotoxicity"; the match must start at a word boundary
DEFAULT_TOPICS: Dict[str, List[str]] = {
    "Safety & toxicity": ["adverse event", "toxicity", "hepatotox", "safety", "mutagen",
                          "carcinogen", "genotox"],
    "Efficacy & clinical trials": ["efficacy", "randomi", "clinical trial", "phase 2",
                                   "phase 3", "placebo", "endpoint"],
    "Pharmacokinetics": ["pharmacokinetic", "bioavailab", "clearance", "half-life",
                         "drug-drug interaction", "exposure"],
    "AI / machine learning": ["machine learning", "deep learning", "neural network",
                              "artificial intelligence", "random forest", "language model"],
    "Biomarkers": ["biomarker", "signature", "prognostic", "predictive marker"],
    "Mechanism & pathways": ["mechanism", "pathway", "signalling", "signaling",
                             "transcription factor", "receptor"],
}


def load_topics(path: str | None) -> Dict[str, List[str]]:
    """Topics from a JSON file ({"Topic name": ["keyword", ...]}) or the built-in defaults."""
    if not path:
        return DEFAULT_TOPICS
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


# --------------------------------------------------------------------------- dedupe
def _norm_title(title: str) -> str:
    return re.sub(r"[^a-z0-9]", "", title.lower())


def _key(rec: dict) -> str:
    if rec.get("doi"):
        return "doi:" + rec["doi"].lower()
    if rec.get("pmid"):
        return "pmid:" + str(rec["pmid"])
    return "title:" + _norm_title(rec.get("title", ""))


def _richness(rec: dict) -> Tuple[int, int]:
    return (1 if rec.get("abstract") else 0, len(rec.get("abstract", "")))


def dedupe(records: List[dict]) -> Tuple[List[dict], int]:
    """Keep one record per paper (the one with the fuller abstract). Returns (records, n_removed)."""
    best: Dict[str, dict] = {}
    for rec in records:
        k = _key(rec)
        if k not in best:
            best[k] = dict(rec, also_pmids=[])
        else:
            kept = best[k]
            if _richness(rec) > _richness(kept):
                also = kept["also_pmids"] + [kept["pmid"]]
                best[k] = dict(rec, also_pmids=also)
            else:
                kept["also_pmids"].append(rec["pmid"])
    unique = list(best.values())
    return unique, len(records) - len(unique)


# --------------------------------------------------------------------------- tagging
def _compile(topics: Dict[str, List[str]]) -> Dict[str, re.Pattern]:
    return {name: re.compile(r"\b(?:" + "|".join(re.escape(k) for k in kws) + ")", re.I)
            for name, kws in topics.items()}


def tag_topics(rec: dict, patterns: Dict[str, re.Pattern]) -> List[str]:
    haystack = " ".join([rec.get("title", ""), rec.get("abstract", ""),
                         " ".join(rec.get("mesh", []))])
    return [name for name, pat in patterns.items() if pat.search(haystack)]


def to_dataframe(records: List[dict], topics: Dict[str, List[str]] = DEFAULT_TOPICS) -> pd.DataFrame:
    patterns = _compile(topics)
    df = pd.DataFrame(records)
    for col in ["pmid", "title", "abstract", "journal", "doi", "url"]:
        if col not in df:
            df[col] = ""
    df["topics"] = [tag_topics(r, patterns) for r in records]
    df["year"] = pd.to_numeric(df.get("year"), errors="coerce").astype("Int64")
    return df.sort_values(["year", "title"], ascending=[False, True]).reset_index(drop=True)


# --------------------------------------------------------------------------- summary
def summarize(df: pd.DataFrame, topics: Dict[str, List[str]] = DEFAULT_TOPICS) -> dict:
    n = len(df)
    years = df["year"].dropna().astype(int)
    year_counts = years.value_counts().sort_index()

    topic_counts = pd.Series({t: int(df["topics"].apply(lambda x, t=t: t in x).sum())
                              for t in topics}).sort_values(ascending=False)

    rows = [{"year": int(y), "topic": t} for y, ts in zip(df["year"], df["topics"])
            if pd.notna(y) for t in ts]
    topic_year = (pd.DataFrame(rows).groupby(["topic", "year"]).size().unstack(fill_value=0)
                  if rows else pd.DataFrame())

    mesh = pd.Series([m for ms in df.get("mesh", []) if isinstance(ms, list) for m in ms])
    journals = df.loc[df["journal"] != "", "journal"].value_counts()

    return {
        "n_records": n,
        "year_range": (int(years.min()), int(years.max())) if len(years) else None,
        "year_counts": year_counts,
        "top_journals": journals.head(10),
        "topic_counts": topic_counts,
        "topic_year": topic_year,
        "top_mesh": mesh.value_counts().head(15) if len(mesh) else mesh,
        "pct_abstract": round(100 * (df["abstract"] != "").mean(), 1) if n else 0.0,
        "pct_doi": round(100 * (df["doi"] != "").mean(), 1) if n else 0.0,
        "n_untagged": int(df["topics"].apply(len).eq(0).sum()),
    }
