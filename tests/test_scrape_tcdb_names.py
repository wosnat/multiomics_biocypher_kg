"""Unit tests for the TCDB names scraper (T6).

Spec: docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md
Plan: docs/superpowers/plans/2026-08-18-tcdb-node-names.md
"""
from multiomics_kg.download import scrape_tcdb_names as m


class TestCleanHtmlFragment:
    def test_strips_tags_and_entities(self):
        s = "<p>Drug:H<sup>+</sup> Antiporter &alpha; &nbsp;porin</p>"
        assert m.clean_html_fragment(s) == "Drug:H + Antiporter α porin"

    def test_tolerates_dangling_unclosed_tag_at_end(self):
        # observed live: one cell ends with an unclosed "<br /"
        s = "Glycerol uptake permease Stl1.<br /"
        assert m.clean_html_fragment(s) == "Glycerol uptake permease Stl1."

    def test_collapses_whitespace(self):
        assert m.clean_html_fragment("a\n  b\t c") == "a b c"


class TestStripCitations:
    def test_removes_year_parenthetical(self):
        s = "Galactose:H+ symporter (Henderson and Giddens 1977). More text."
        assert m.strip_citations(s) == "Galactose:H+ symporter. More text."

    def test_keeps_non_citation_parenthetical(self):
        s = "The Sugar Porter (SP) Family"
        assert m.strip_citations(s) == "The Sugar Porter (SP) Family"

    def test_removes_multi_ref_group(self):
        s = "GalP (Henderson 1977; Patching et al., 2008). End."
        assert m.strip_citations(s) == "GalP. End."

    def test_nested_citation_removed_content_parenthetical_kept(self):
        # inner group is the citation; the outer group is annotation content
        s = "Stl1 (similar to Stl1 of S. cerevisiae (Ferreira et al., 2005))."
        assert m.strip_citations(s) == "Stl1 (similar to Stl1 of S. cerevisiae)."


class TestFirstSentence:
    def test_plain_split(self):
        s = "MDR efflux pump, MdeA. Exports quaternary ammonium compounds."
        assert m.first_sentence(s) == "MDR efflux pump, MdeA."

    def test_does_not_split_on_genus_initial(self):
        s = "Transporter of E. coli origin. Second sentence."
        assert m.first_sentence(s) == "Transporter of E. coli origin."

    def test_does_not_split_on_et_al(self):
        s = "Characterized by Smith et al. Second clause continues here."
        assert m.first_sentence(s) == "Characterized by Smith et al. Second clause continues here."

    def test_no_terminal_period(self):
        s = "The glucose transport protein, GTP1"
        assert m.first_sentence(s) == "The glucose transport protein, GTP1"


class TestCaps:
    def test_short_string_unchanged(self):
        assert m.cap_at_word_boundary("short name", 150) == "short name"

    def test_long_string_capped_at_word_boundary(self):
        s = ("word " * 60).strip()  # 299 chars
        out = m.cap_at_word_boundary(s, 150)
        assert len(out) <= 151  # 150 + ellipsis char
        assert out.endswith("…")
        assert not out[:-1].endswith(" ")

    def test_make_name_pipeline(self):
        s = ("Glycerol-P:Pi antiporter (Ambudkar et al. 1986). "
             "The 3-d structure is known.")
        assert m.make_name(s) == "Glycerol-P:Pi antiporter."
