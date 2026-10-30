"""Offline tests: run with  pytest -q  (no network needed)."""
import json
import os

import pytest

from analyze_literature import DEFAULT_TOPICS, dedupe, summarize, tag_topics, to_dataframe, _compile
from build_dashboard import build, pdf_safe
from pubmed_client import fetch_records, parse_pubmed_xml, search_pmids

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = os.path.join(HERE, "sample_data", "pubmed_sample.xml")


@pytest.fixture(scope="module")
def records():
    with open(SAMPLE, encoding="utf-8") as fh:
        return parse_pubmed_xml(fh.read())


# ---------------------------------------------------------------- parsing
def test_parses_all_articles(records):
    assert len(records) == 63            # 60 papers + 3 duplicate-DOI records


def test_title_inline_tags_are_flattened(records):
    assert all("<" not in r["title"] for r in records)
    assert any("In vitro" in r["title"] for r in records)


def test_special_characters_survive(records):
    assert any("TGF-β" in r["title"] and "&" in r["title"] for r in records)


def test_structured_abstract_labels(records):
    structured = [r for r in records if r["abstract"].startswith("Background:")]
    assert structured and "Methods:" in structured[0]["abstract"]


def test_missing_abstract_and_doi_are_empty_strings(records):
    assert any(r["abstract"] == "" for r in records)
    assert any(r["doi"] == "" for r in records)


def test_medlinedate_year_is_parsed(records):
    assert all(isinstance(r["year"], int) and 2018 <= r["year"] <= 2025 for r in records)


# ---------------------------------------------------------------- dedupe / tagging / summary
def test_dedupe_removes_duplicates_and_keeps_fuller_record(records):
    unique, removed = dedupe(records)
    assert removed == 3 and len(unique) == 60
    assert all(r["abstract"] or r["pmid"] not in {"90000061", "90000062", "90000063"} for r in unique)
    assert any(r["also_pmids"] for r in unique)


def test_topic_tagging_prefix_match():
    pats = _compile(DEFAULT_TOPICS)
    rec = {"title": "Hepatotoxicity in trial", "abstract": "", "mesh": []}
    assert "Safety & toxicity" in tag_topics(rec, pats)
    assert "AI / machine learning" not in tag_topics(rec, pats)


def test_topic_requires_word_start():
    pats = _compile({"T": ["safety"]})
    assert tag_topics({"title": "unsafety", "abstract": "", "mesh": []}, pats) == []


def test_summary_is_consistent(records):
    unique, _ = dedupe(records)
    s = summarize(to_dataframe(unique))
    assert s["n_records"] == 60
    assert int(s["year_counts"].sum()) == 60
    assert 0 <= s["pct_abstract"] <= 100


# ---------------------------------------------------------------- client with a fake HTTP layer
def test_client_batches_requests_and_parses():
    calls = []

    def fake(url, params):
        calls.append(url.rsplit("/", 1)[-1])
        if "esearch" in url:
            return json.dumps({"esearchresult": {"idlist": [str(1000 + i) for i in range(450)]}})
        ids = params["id"].split(",")
        arts = "".join(f"<PubmedArticle><MedlineCitation><PMID>{i}</PMID><Article>"
                       f"<ArticleTitle>T{i}</ArticleTitle></Article></MedlineCitation></PubmedArticle>"
                       for i in ids)
        return f"<PubmedArticleSet>{arts}</PubmedArticleSet>"

    pmids = search_pmids("anything", retmax=450, fetch=fake)
    recs = fetch_records(pmids, batch_size=200, delay=0, fetch=fake)
    assert len(pmids) == 450 and len(recs) == 450
    assert calls == ["esearch.fcgi"] + ["efetch.fcgi"] * 3          # 200 + 200 + 50


def test_search_adds_date_range():
    seen = {}

    def fake(url, params):
        seen.update(params)
        return json.dumps({"esearchresult": {"idlist": []}})

    search_pmids("x", mindate="2020", fetch=fake)
    assert seen["datetype"] == "pdat" and seen["mindate"] == "2020"


# ---------------------------------------------------------------- outputs
def test_pdf_safe_maps_symbols():
    assert pdf_safe("TGF-β ≥ 2") == "TGF-beta >= 2"


def test_build_creates_dashboard_and_pdf(tmp_path):
    summary = build(SAMPLE, str(tmp_path), "Test title")
    html = (tmp_path / "dashboard.html").read_text(encoding="utf-8")
    pdf = (tmp_path / "literature_report.pdf").read_bytes()
    assert summary["n_records"] == 60
    assert "Test title" in html and "const DATA" in html
    assert html.count("</script>") == 1            # data cannot break out of the script tag
    assert pdf.startswith(b"%PDF") and len(pdf) > 10_000
