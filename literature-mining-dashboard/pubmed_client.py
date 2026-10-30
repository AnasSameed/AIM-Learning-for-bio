"""
PubMed client built on NCBI E-utilities (esearch + efetch).

    python pubmed_client.py "nitrosamine AND mutagenicity" --max 200 --out data/records.json

What it does
  1. esearch  -> list of PMIDs for a query
  2. efetch   -> PubMed XML for those PMIDs, fetched in batches
  3. parse    -> clean records: pmid, title, abstract, journal, year, doi, authors, mesh, pub_types

Design notes
  * The HTTP layer is injectable (`fetch=`), so the whole pipeline can be tested offline.
  * NCBI allows ~3 requests/second without an API key (10/s with one). The client sleeps
    between batches and reads NCBI_API_KEY from the environment if it is set.
  * esearch returns at most 10,000 ids per query; narrow with --mindate/--maxdate for more.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import Callable, Dict, Iterable, List, Optional

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "aim-literature-dashboard"

Fetcher = Callable[[str, Dict[str, str]], str]


# --------------------------------------------------------------------------- HTTP
def http_get(url: str, params: Dict[str, str], retries: int = 4, timeout: int = 30) -> str:
    """GET with simple exponential backoff on rate-limit / server errors."""
    import requests  # imported here so parsing/tests work without requests installed

    for attempt in range(retries):
        resp = requests.get(url, params=params, timeout=timeout)
        if resp.status_code == 200:
            return resp.text
        if resp.status_code in (429, 500, 502, 503, 504):
            time.sleep(2 ** attempt)
            continue
        resp.raise_for_status()
    raise RuntimeError(f"NCBI request failed after {retries} attempts: {url}")


def _base_params(email: Optional[str], api_key: Optional[str]) -> Dict[str, str]:
    params = {"tool": TOOL_NAME}
    email = email or os.environ.get("NCBI_EMAIL")
    api_key = api_key or os.environ.get("NCBI_API_KEY")
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    return params


# --------------------------------------------------------------------------- search
def search_pmids(
    term: str,
    retmax: int = 200,
    mindate: Optional[str] = None,
    maxdate: Optional[str] = None,
    fetch: Fetcher = http_get,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> List[str]:
    """Return PMIDs (as strings) matching a PubMed query, newest first."""
    params = _base_params(email, api_key)
    params.update({"db": "pubmed", "term": term, "retmode": "json",
                   "retmax": str(retmax), "sort": "pub_date"})
    if mindate or maxdate:
        params.update({"datetype": "pdat",
                       "mindate": mindate or "1900",
                       "maxdate": maxdate or "3000"})
    payload = json.loads(fetch(f"{EUTILS}/esearch.fcgi", params))
    return [str(i) for i in payload.get("esearchresult", {}).get("idlist", [])]


# --------------------------------------------------------------------------- fetch
def _chunks(items: List[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield items[i:i + size]


def fetch_records(
    pmids: List[str],
    batch_size: int = 200,
    delay: float = 0.34,
    fetch: Fetcher = http_get,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> List[dict]:
    """Download and parse PubMed records for a list of PMIDs, in batches."""
    records: List[dict] = []
    for n, batch in enumerate(_chunks(list(pmids), batch_size)):
        if n and delay:
            time.sleep(delay)
        params = _base_params(email, api_key)
        params.update({"db": "pubmed", "id": ",".join(batch),
                       "retmode": "xml", "rettype": "abstract"})
        records.extend(parse_pubmed_xml(fetch(f"{EUTILS}/efetch.fcgi", params)))
    return records


# --------------------------------------------------------------------------- parsing
_WS = re.compile(r"\s+")
_YEAR = re.compile(r"(?:19|20)\d{2}")


def _text(el: Optional[ET.Element]) -> str:
    """All text inside an element (keeps text inside inline tags such as <i>, <sub>)."""
    if el is None:
        return ""
    return _WS.sub(" ", "".join(el.itertext())).strip()


def _year(article: ET.Element) -> Optional[int]:
    pub = article.find("./MedlineCitation/Article/Journal/JournalIssue/PubDate")
    if pub is not None:
        y = pub.findtext("Year")
        if y and y.isdigit():
            return int(y)
        m = _YEAR.search(pub.findtext("MedlineDate") or "")
        if m:
            return int(m.group())
    y = article.findtext("./MedlineCitation/Article/ArticleDate/Year")
    return int(y) if y and y.isdigit() else None


def _abstract(article: ET.Element) -> str:
    parts = []
    for node in article.findall("./MedlineCitation/Article/Abstract/AbstractText"):
        label = node.get("Label")
        body = _text(node)
        if body:
            parts.append(f"{label.title()}: {body}" if label else body)
    return " ".join(parts)


def _doi(article: ET.Element) -> str:
    for aid in article.findall("./PubmedData/ArticleIdList/ArticleId"):
        if aid.get("IdType") == "doi" and aid.text:
            return aid.text.strip().lower()
    for loc in article.findall("./MedlineCitation/Article/ELocationID"):
        if loc.get("EIdType") == "doi" and loc.text:
            return loc.text.strip().lower()
    return ""


def _authors(article: ET.Element) -> List[str]:
    names = []
    for a in article.findall("./MedlineCitation/Article/AuthorList/Author"):
        collective = a.findtext("CollectiveName")
        if collective:
            names.append(collective.strip())
            continue
        last, initials = a.findtext("LastName"), a.findtext("Initials")
        if last:
            names.append(f"{last} {initials}".strip() if initials else last)
    return names


def parse_pubmed_xml(xml_text: str) -> List[dict]:
    """Turn a PubMed efetch XML document into a list of flat record dicts."""
    root = ET.fromstring(xml_text)
    records = []
    for art in root.iter("PubmedArticle"):
        pmid = (art.findtext("./MedlineCitation/PMID") or "").strip()
        if not pmid:
            continue
        records.append({
            "pmid": pmid,
            "title": _text(art.find("./MedlineCitation/Article/ArticleTitle")).rstrip("."),
            "abstract": _abstract(art),
            "journal": _text(art.find("./MedlineCitation/Article/Journal/Title")),
            "year": _year(art),
            "doi": _doi(art),
            "authors": _authors(art),
            "mesh": [_text(d) for d in art.findall(
                "./MedlineCitation/MeshHeadingList/MeshHeading/DescriptorName")],
            "pub_types": [_text(p) for p in art.findall(
                "./MedlineCitation/Article/PublicationTypeList/PublicationType")],
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })
    return records


# --------------------------------------------------------------------------- CLI
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Search PubMed and save parsed records as JSON.")
    ap.add_argument("term", help='PubMed query, e.g. "nitrosamine AND mutagenicity"')
    ap.add_argument("--max", type=int, default=200, help="maximum number of records (default 200)")
    ap.add_argument("--mindate", help="earliest publication date, e.g. 2020 or 2020/06")
    ap.add_argument("--maxdate", help="latest publication date")
    ap.add_argument("--email", help="contact email sent to NCBI (or set NCBI_EMAIL)")
    ap.add_argument("--out", default="data/records.json", help="output JSON file")
    args = ap.parse_args(argv)

    pmids = search_pmids(args.term, args.max, args.mindate, args.maxdate, email=args.email)
    print(f"{len(pmids)} PMIDs found")
    records = fetch_records(pmids, email=args.email)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(records, fh, indent=2, ensure_ascii=False)
    print(f"{len(records)} records saved to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
