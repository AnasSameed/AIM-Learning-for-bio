"""
Build a self-contained HTML dashboard and a PDF reference report from literature records.

    python build_dashboard.py --input sample_data/pubmed_sample.xml --out-dir output
    python build_dashboard.py --input data/records.json --topics topics.json --title "Nitrosamine literature"

Inputs : PubMed efetch XML, or the JSON written by pubmed_client.py
Outputs: <out-dir>/dashboard.html          searchable, filterable, one file, no server needed
         <out-dir>/literature_report.pdf   summary charts + reference digest grouped by topic
"""
from __future__ import annotations

import argparse
import base64
import datetime as dt
import io
import json
import os
import sys
import unicodedata
from typing import Dict, List
from xml.sax.saxutils import escape

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from jinja2 import Template  # noqa: E402

from analyze_literature import dedupe, load_topics, summarize, to_dataframe  # noqa: E402
from pubmed_client import parse_pubmed_xml  # noqa: E402

BLUE = "#1D6FA5"


# --------------------------------------------------------------------------- charts
def _png(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def _clean(ax, axis="y"):
    """Hide top/right spines and force whole-number ticks on the count axis."""
    from matplotlib.ticker import MaxNLocator

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_locator(MaxNLocator(integer=True))


def make_charts(summary: dict) -> Dict[str, bytes]:
    charts: Dict[str, bytes] = {}

    yc = summary["year_counts"]
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.bar(yc.index.astype(str), yc.values, color=BLUE)
    ax.set_title("Publications per year", loc="left", fontsize=11)
    ax.set_ylabel("Papers")
    _clean(ax)
    charts["years"] = _png(fig)

    tc = summary["topic_counts"].sort_values()
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.barh(tc.index, tc.values, color=BLUE)
    ax.set_title("Papers per topic (a paper can have several)", loc="left", fontsize=11)
    _clean(ax, axis="x")
    charts["topics"] = _png(fig)

    tj = summary["top_journals"].sort_values()
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.barh([j if len(j) <= 38 else j[:36] + ".." for j in tj.index], tj.values, color=BLUE)
    ax.set_title("Top journals", loc="left", fontsize=11)
    _clean(ax, axis="x")
    charts["journals"] = _png(fig)

    ty = summary["topic_year"]
    if not ty.empty:
        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        ax.imshow(ty.values, cmap="Blues", aspect="auto")
        ax.set_xticks(range(len(ty.columns)), [str(c) for c in ty.columns])
        ax.set_yticks(range(len(ty.index)), ty.index)
        for i in range(ty.shape[0]):
            for j in range(ty.shape[1]):
                v = ty.values[i, j]
                ax.text(j, i, v, ha="center", va="center", fontsize=8,
                        color="white" if v > ty.values.max() * 0.6 else "black")
        ax.set_title("Topic trend by year", loc="left", fontsize=11)
        charts["heatmap"] = _png(fig)
    return charts


# --------------------------------------------------------------------------- HTML
HTML_TEMPLATE = Template("""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{{ title }}</title>
<style>
 body{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:0;background:#f5f7fa;color:#1c2733}
 header{background:#12395a;color:#fff;padding:22px 28px}
 header h1{margin:0;font-size:22px} header p{margin:4px 0 0;opacity:.8;font-size:13px}
 main{max-width:1100px;margin:0 auto;padding:20px}
 .kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-bottom:18px}
 .kpi{background:#fff;border-radius:10px;padding:14px;box-shadow:0 1px 3px rgba(0,0,0,.08)}
 .kpi b{display:block;font-size:24px;color:#12395a} .kpi span{font-size:12px;color:#5b6b7b}
 .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:14px}
 .card{background:#fff;border-radius:10px;padding:10px;box-shadow:0 1px 3px rgba(0,0,0,.08)}
 .card img{width:100%;height:auto}
 .controls{display:flex;gap:10px;flex-wrap:wrap;margin:22px 0 10px}
 input,select{padding:8px 10px;border:1px solid #c5cfd9;border-radius:6px;font-size:14px}
 input{flex:1;min-width:220px}
 table{width:100%;border-collapse:collapse;background:#fff;border-radius:10px;overflow:hidden;font-size:14px}
 th,td{text-align:left;padding:9px 10px;border-bottom:1px solid #e6ebf0;vertical-align:top}
 th{background:#eaf0f6;cursor:pointer;white-space:nowrap}
 .tag{display:inline-block;background:#e3eef8;color:#12395a;border-radius:10px;padding:1px 8px;margin:1px 3px 1px 0;font-size:11px}
 details summary{cursor:pointer;color:#1d6fa5;font-size:12px} details p{margin:6px 0 0;font-size:13px;color:#3a4856}
 footer{color:#6b7a89;font-size:12px;padding:20px 0}
</style></head><body>
<header><h1>{{ title }}</h1><p>Generated {{ generated }} &middot; {{ n_removed }} duplicate record(s) removed</p></header>
<main>
 <div class="kpis">
  <div class="kpi"><b>{{ s.n_records }}</b><span>unique papers</span></div>
  <div class="kpi"><b>{{ years }}</b><span>year range</span></div>
  <div class="kpi"><b>{{ s.pct_abstract }}%</b><span>with abstract</span></div>
  <div class="kpi"><b>{{ s.pct_doi }}%</b><span>with DOI</span></div>
  <div class="kpi"><b>{{ s.n_untagged }}</b><span>not matched to a topic</span></div>
 </div>
 <div class="grid">{% for key, img in charts.items() %}<div class="card"><img alt="{{ key }} chart" src="data:image/png;base64,{{ img }}"></div>{% endfor %}</div>
 <div class="controls">
  <input id="q" placeholder="Search title, abstract, authors, journal...">
  <select id="topic"><option value="">All topics</option>{% for t in topics %}<option>{{ t }}</option>{% endfor %}</select>
 </div>
 <table><thead><tr><th data-k="year">Year</th><th data-k="title">Paper</th><th data-k="journal">Journal</th><th>Topics</th></tr></thead><tbody id="rows"></tbody></table>
 <p id="count"></p>
 <footer>Topic tags come from keyword matching and are a first-pass screen; review before relying on them.</footer>
</main>
<script>
const DATA = {{ data_json }};
let sortKey = "year", sortDir = -1;
const $ = id => document.getElementById(id);
function render(){
  const q = $("q").value.toLowerCase(), t = $("topic").value;
  const rows = DATA.filter(r => (!t || r.topics.includes(t)) &&
      (!q || (r.title+" "+r.abstract+" "+r.authors+" "+r.journal).toLowerCase().includes(q)))
    .sort((a,b) => (a[sortKey] > b[sortKey] ? 1 : a[sortKey] < b[sortKey] ? -1 : 0) * sortDir);
  const body = $("rows"); body.textContent = "";
  rows.forEach(r => {
    const tr = document.createElement("tr");
    const y = document.createElement("td"); y.textContent = r.year || ""; tr.appendChild(y);
    const p = document.createElement("td");
    const a = document.createElement("a"); a.href = r.url; a.target = "_blank"; a.rel = "noopener"; a.textContent = r.title;
    p.appendChild(a);
    const meta = document.createElement("div"); meta.style.cssText = "font-size:12px;color:#5b6b7b";
    meta.textContent = r.authors + (r.doi ? "  |  doi:" + r.doi : "") + "  |  PMID " + r.pmid; p.appendChild(meta);
    if (r.abstract){ const d = document.createElement("details"), s = document.createElement("summary"), x = document.createElement("p");
      s.textContent = "Abstract"; x.textContent = r.abstract; d.appendChild(s); d.appendChild(x); p.appendChild(d); }
    tr.appendChild(p);
    const j = document.createElement("td"); j.textContent = r.journal; tr.appendChild(j);
    const tg = document.createElement("td");
    r.topics.forEach(x => { const s = document.createElement("span"); s.className = "tag"; s.textContent = x; tg.appendChild(s); });
    tr.appendChild(tg); body.appendChild(tr);
  });
  $("count").textContent = rows.length + " of " + DATA.length + " papers shown";
}
$("q").addEventListener("input", render); $("topic").addEventListener("change", render);
document.querySelectorAll("th[data-k]").forEach(th => th.addEventListener("click", () => {
  const k = th.dataset.k; sortDir = (sortKey === k) ? -sortDir : 1; sortKey = k; render(); }));
render();
</script></body></html>""")


def _authors_short(authors: List[str]) -> str:
    return ", ".join(authors[:3]) + (" et al." if len(authors) > 3 else "")


def render_html(df, summary: dict, charts: Dict[str, bytes], title: str, out_path: str) -> None:
    data = [{"pmid": r.pmid, "title": r.title, "abstract": r.abstract, "journal": r.journal,
             "year": int(r.year) if r.year is not None and str(r.year) != "<NA>" else 0,
             "doi": r.doi, "url": r.url, "authors": _authors_short(list(r.authors or [])),
             "topics": list(r.topics)} for r in df.itertuples()]
    html = HTML_TEMPLATE.render(
        title=title, generated=dt.date.today().isoformat(), s=summary,
        n_removed=summary.get("n_removed", 0),
        years="-".join(map(str, summary["year_range"])) if summary["year_range"] else "n/a",
        topics=list(summary["topic_counts"].index),
        charts={k: base64.b64encode(v).decode() for k, v in charts.items()},
        # "</" is escaped so abstract text can never close the script tag
        data_json=json.dumps(data, ensure_ascii=False).replace("</", "<\\/"))
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)


# --------------------------------------------------------------------------- PDF
_GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "κ": "kappa",
          "μ": "u", "≥": ">=", "≤": "<=", "→": "->", "−": "-"}


def pdf_safe(text: str) -> str:
    """The built-in PDF fonts only cover Latin-1/cp1252; map common symbols and drop the rest."""
    text = unicodedata.normalize("NFC", text or "")
    for src, dst in _GREEK.items():
        text = text.replace(src, dst)
    return text.encode("cp1252", "replace").decode("cp1252")


def render_pdf(df, summary: dict, charts: Dict[str, bytes], title: str, out_path: str) -> None:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph,
                                    SimpleDocTemplate, Spacer, Table, TableStyle)

    ss = getSampleStyleSheet()
    h1 = ParagraphStyle("h1", parent=ss["Title"], fontSize=20, leading=24, alignment=0)
    h2 = ParagraphStyle("h2", parent=ss["Heading2"], textColor=colors.HexColor("#12395a"), spaceBefore=10)
    body = ParagraphStyle("b", parent=ss["BodyText"], fontSize=9.5, leading=12.5)
    small = ParagraphStyle("s", parent=body, fontSize=8.5, leading=11, textColor=colors.HexColor("#444444"))
    esc = lambda s: escape(pdf_safe(s))  # noqa: E731

    story = [Paragraph(esc(title), h1),
             Paragraph(f"Literature report &middot; generated {dt.date.today().isoformat()}", small),
             Spacer(1, 10)]
    yr = "-".join(map(str, summary["year_range"])) if summary["year_range"] else "n/a"
    kpi = Table([[str(summary["n_records"]), yr, f"{summary['pct_abstract']}%", f"{summary['pct_doi']}%"],
                 ["unique papers", "year range", "with abstract", "with DOI"]], colWidths=[4.2 * cm] * 4)
    kpi.setStyle(TableStyle([("FONTSIZE", (0, 0), (-1, 0), 15), ("FONTSIZE", (0, 1), (-1, 1), 8),
                             ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#12395a")),
                             ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#eaf0f6")),
                             ("ALIGN", (0, 0), (-1, -1), "CENTER")]))
    story += [kpi, Spacer(1, 10)]
    for key in ("years", "topics", "journals", "heatmap"):
        if key in charts:
            story.append(Image(io.BytesIO(charts[key]), width=15.5 * cm, height=7.6 * cm))
            story.append(Spacer(1, 6))
    story.append(PageBreak())

    story.append(Paragraph("Reference digest (grouped by first matching topic)", h2))
    story.append(Paragraph("Topic tags are keyword matches, a first-pass screen to be reviewed by a person.", small))
    groups: Dict[str, list] = {}
    for r in df.itertuples():
        groups.setdefault(r.topics[0] if r.topics else "Not matched to a topic", []).append(r)
    for topic, items in groups.items():
        story.append(Paragraph(f"{esc(topic)} ({len(items)})", h2))
        for r in items:
            abstract = (r.abstract[:330] + "...") if len(r.abstract) > 330 else r.abstract
            ids = f"PMID {r.pmid}" + (f" | doi:{r.doi}" if r.doi else "")
            entry = [Paragraph(f"<b>{esc(r.title)}</b>", body),
                     Paragraph(esc(f"{_authors_short(list(r.authors or [])).rstrip('.')}. {r.journal}, "
                                   f"{r.year if str(r.year) != '<NA>' else 'n.d.'}. {ids}"), small)]
            if abstract:
                entry.append(Paragraph(esc(abstract), small))
            entry.append(Spacer(1, 5))
            story.append(KeepTogether(entry))          # never split a reference across pages

    def footer(canvas, doc):
        canvas.setFont("Helvetica", 8)
        canvas.drawString(2 * cm, 1.2 * cm, pdf_safe(title)[:80])
        canvas.drawRightString(A4[0] - 2 * cm, 1.2 * cm, f"Page {doc.page}")

    SimpleDocTemplate(out_path, pagesize=A4, leftMargin=2 * cm, rightMargin=2 * cm,
                      topMargin=1.8 * cm, bottomMargin=2 * cm, title=pdf_safe(title)
                      ).build(story, onFirstPage=footer, onLaterPages=footer)


# --------------------------------------------------------------------------- CLI
def load_records(path: str) -> List[dict]:
    if path.lower().endswith(".xml"):
        with open(path, encoding="utf-8") as fh:
            return parse_pubmed_xml(fh.read())
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def build(input_path: str, out_dir: str, title: str, topics_path: str | None = None) -> dict:
    topics = load_topics(topics_path)
    records, n_removed = dedupe(load_records(input_path))
    df = to_dataframe(records, topics)
    summary = summarize(df, topics)
    summary["n_removed"] = n_removed
    charts = make_charts(summary)
    os.makedirs(out_dir, exist_ok=True)
    render_html(df, summary, charts, title, os.path.join(out_dir, "dashboard.html"))
    render_pdf(df, summary, charts, title, os.path.join(out_dir, "literature_report.pdf"))
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the literature dashboard (HTML) and report (PDF).")
    ap.add_argument("--input", required=True, help="PubMed XML or records.json")
    ap.add_argument("--out-dir", default="output")
    ap.add_argument("--topics", help="JSON file of topic keywords (optional)")
    ap.add_argument("--title", default="Literature dashboard")
    args = ap.parse_args(argv)
    s = build(args.input, args.out_dir, args.title, args.topics)
    print(f"{s['n_records']} unique papers ({s['n_removed']} duplicates removed) -> {args.out_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
