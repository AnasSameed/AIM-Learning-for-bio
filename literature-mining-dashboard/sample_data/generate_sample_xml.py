"""
Generate a SYNTHETIC PubMed-format XML file for testing and demos.

Nothing here is a real publication: every title starts with "[SAMPLE]", PMIDs are in the
9xxxxxxx range, DOIs use a made-up prefix and journals are fictional. The structure
(structured abstracts, MedlineDate, inline <i> tags, missing abstracts/DOIs, duplicate DOIs)
mimics the edge cases of real PubMed XML so the parser and de-duplication can be tested.

    python generate_sample_xml.py            # writes pubmed_sample.xml next to this script
"""
import os
import random
from xml.sax.saxutils import escape

random.seed(2026)

JOURNALS = ["Journal of Sample Pharmacology", "Annals of Demo Medicine",
            "Synthetic Biology Letters", "Archive of Mock Toxicology",
            "Example Clinical Reports", "Placeholder Bioinformatics"]

# (topic phrase, keywords that the default topic tagger should pick up)
THEMES = [
    ("hepatotoxicity signal detection", "safety and hepatotoxicity were assessed"),
    ("randomised placebo-controlled trial", "a randomised placebo-controlled phase 3 clinical trial measured efficacy"),
    ("pharmacokinetic modelling", "pharmacokinetic parameters including clearance and half-life were estimated"),
    ("machine learning prediction", "a random forest and deep learning model were trained"),
    ("biomarker discovery", "a prognostic biomarker signature was identified"),
    ("receptor signalling pathway", "the mechanism of receptor signalling was dissected"),
    ("mutagenicity assessment", "genotoxicity and mutagenic potency were compared"),
]
CONDITIONS = ["prostate cancer", "type 2 diabetes", "rheumatoid arthritis",
              "Alzheimer disease", "hypertension", "asthma"]
MESH = ["Humans", "Animals", "Mice", "Biomarkers", "Machine Learning", "Neoplasms",
        "Drug-Related Side Effects and Adverse Reactions", "Pharmacokinetics",
        "Randomized Controlled Trials as Topic", "Signal Transduction"]
LAST = ["Rao", "Mehta", "Iyer", "Khan", "Nair", "Sharma", "Das", "Patel", "Singh", "Verma"]


def article(pmid, title_xml, abstract_xml, journal, year_xml, doi, authors, mesh):
    abs_block = f"<Abstract>{abstract_xml}</Abstract>" if abstract_xml else ""
    doi_block = f'<ArticleId IdType="doi">{doi}</ArticleId>' if doi else ""
    auth = "".join(f"<Author><LastName>{l}</LastName><Initials>{i}</Initials></Author>"
                   for l, i in authors)
    mh = "".join(f"<MeshHeading><DescriptorName>{escape(m)}</DescriptorName></MeshHeading>" for m in mesh)
    return f"""<PubmedArticle><MedlineCitation><PMID>{pmid}</PMID><Article>
<Journal><JournalIssue><PubDate>{year_xml}</PubDate></JournalIssue><Title>{escape(journal)}</Title></Journal>
<ArticleTitle>{title_xml}</ArticleTitle>{abs_block}<AuthorList>{auth}</AuthorList>
<PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList>
</Article><MeshHeadingList>{mh}</MeshHeadingList></MedlineCitation>
<PubmedData><ArticleIdList><ArticleId IdType="pubmed">{pmid}</ArticleId>{doi_block}</ArticleIdList></PubmedData></PubmedArticle>"""


def build(n=60):
    items, pmid = [], 90000001
    for i in range(n):
        theme, sentence = random.choice(THEMES)
        cond = random.choice(CONDITIONS)
        year = random.choices(range(2018, 2026), weights=[1, 1, 2, 3, 4, 5, 6, 7])[0]
        title = f"[SAMPLE] {theme.capitalize()} in {cond}"
        if i % 9 == 0:                                   # inline markup inside a title
            title = f"[SAMPLE] <i>In vitro</i> {theme} in {cond}"
        if i % 11 == 0:                                  # special characters
            title += " &amp; TGF-β signalling"
        body = f"In patients with {cond}, {sentence}."
        if i % 4 == 0:                                   # structured abstract
            abstract = (f'<AbstractText Label="BACKGROUND">{escape(cond.capitalize())} is a major burden.</AbstractText>'
                        f'<AbstractText Label="METHODS">{escape(sentence.capitalize())}.</AbstractText>'
                        f'<AbstractText Label="RESULTS">Synthetic results for demonstration only.</AbstractText>')
        elif i % 13 == 0:                                # no abstract
            abstract = ""
        else:
            abstract = f"<AbstractText>{escape(body)} Synthetic text for demonstration only.</AbstractText>"
        year_xml = f"<Year>{year}</Year>"
        if i % 17 == 0:                                  # MedlineDate instead of Year
            year_xml = f"<MedlineDate>{year} Jan-Feb</MedlineDate>"
        doi = "" if i % 15 == 0 else f"10.0000/sample.{1000 + i}"
        authors = [(random.choice(LAST), random.choice("ABCDKMNRSV")) for _ in range(random.randint(2, 6))]
        mesh = random.sample(MESH, 3)
        items.append((pmid, title, abstract, random.choice(JOURNALS), year_xml, doi, authors, mesh))
        pmid += 1
    # three duplicates: same DOI as an earlier paper but a new PMID (e.g. preprint vs published)
    for src in (1, 2, 3):
        p, t, a, j, y, d, au, m = items[src]
        items.append((pmid, t, "", j, y, d, au, m))      # duplicate without abstract
        pmid += 1
    body = "".join(article(*it) for it in items)
    return f'<?xml version="1.0" encoding="UTF-8"?>\n<PubmedArticleSet>{body}</PubmedArticleSet>\n'


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pubmed_sample.xml")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(build())
    print("wrote", out)
