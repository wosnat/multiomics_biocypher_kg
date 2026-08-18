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


BROWSE_FIXTURE = """
<ul id="red" class="treeview-red">
  <li><span>
    <div rel="1" class="entry">
      <div class="tcid name">&nbsp;1:&nbsp;Channels/Pores</div>
    </div>
  </span>
  <ul><li><span>
    <div rel="1.A" class="entry">
      <div class="tcid name">&nbsp;1.A:&nbsp;&#945;-Type Channels</div>
    </div>
  </span>
  <ul><li><span>
    <div rel="1.A.1" class="entry" style="border-bottom:1px dotted #999;">
      <div class="tcid name" style="cursor:pointer;">&nbsp;1.A.1:&nbsp;The Voltage-gated Ion Channel (VIC) Superfamily </div>
    </div>
  </span></li></ul>
  </li></ul></li>
</ul>
"""

FAMILY_PAGE_FIXTURE = """
<A HREF="/search/result.php?tc=2.A.1">Read Family Description</A>
<table id="result-cluster" style="vertical-align:top">
  <tr><th>TCID</th><th>Name</th><th>Domain</th><th>Kingdom/Phylum</th><th class="right-border">Protein(s)</th></tr>
  <tr><td colspan="4" id="right-border"><strong><A id="2.A.1.1"></A>2.A.1.1:&nbsp;&nbsp;The Sugar Porter (SP) Family</strong></td></tr>
  <tr><td width="60" style="vertical-align:top">
    <A id="2.A.1.1.1"></A>
    <A valign="top" HREF="/search/result.php?tc=2.A.1.1.1">
2.A.1.1.1</A><br></td><td><div class='400aa753232406c9e9e7370875b4f60b'><p>Galactose:H<sup>+</sup> symporter, GalP (<a class="reflink" href="/search/result.php?tc=2.A.1#ref7167">Hern&aacute;ndez-Montalvo <em>et al.</em>, 2001</a>). Long tail text.</p></div></td>
    <td>Bacteria</td><td>Pseudomonadota</td>
    <td id="right-border" width="170"><div class='x'>GalP of <em>E. coli</em></div></td></tr>
  <tr><td colspan="4" id="right-border"><strong><A id="2.A.1.2"></A>2.A.1.2:&nbsp;&nbsp;The Drug:H<sup>+</sup> Antiporter-1 (12 Spanner) (DHA1) Family</strong></td></tr>
</table>
<div id="sidebar"><a HREF="/search/result.php?tc=9.Z.9.9.9">9.Z.9.9.9</a></div>
"""

UNNAMED_FAMILY_PAGE_FIXTURE = """
<table id="result-cluster" style="vertical-align:top">
  <tr><th>TCID</th><th>Name</th><th>Domain</th><th>Kingdom/Phylum</th><th class="right-border">Protein(s)</th></tr>
  <tr><td width="60" style="vertical-align:top">
    <A id="1.A.11.1"></A>
    <A id="1.A.11.1.1"></A>
    <A valign="top" HREF="/search/result.php?tc=1.A.11.1.1">
1.A.11.1.1</A><br></td><td><div class='cf8b031f73d53100bb24514ad08c11f0'><p>Ammonia transporter and regulatory sensor, AmtB.</p></div></td>
    <td>Bacteria</td><td>Pseudomonadota</td>
    <td id="right-border"><div class='y'>AmtB of <em>E. coli</em></div></td></tr>
</table>
"""


class TestParseBrowse:
    def test_parses_all_three_levels(self):
        names = m.parse_browse(BROWSE_FIXTURE)
        assert names["1"] == "Channels/Pores"
        assert names["1.A"] == "α-Type Channels"
        assert names["1.A.1"] == "The Voltage-gated Ion Channel (VIC) Superfamily"

    def test_no_4part_entries(self):
        assert all(k.count(".") <= 2 for k in m.parse_browse(BROWSE_FIXTURE))


class TestParseFamilyPage:
    def test_named_subfamily_headers(self):
        subfams, _systems = m.parse_family_page(FAMILY_PAGE_FIXTURE)
        assert subfams["2.A.1.1"] == "The Sugar Porter (SP) Family"
        # inline markup inside the name must survive the capture (spec 6a):
        assert subfams["2.A.1.2"] == "The Drug:H + Antiporter-1 (12 Spanner) (DHA1) Family"

    def test_system_full_text(self):
        _, systems = m.parse_family_page(FAMILY_PAGE_FIXTURE)
        assert systems["2.A.1.1.1"].startswith("Galactose:H + symporter, GalP")
        assert "Long tail text." in systems["2.A.1.1.1"]
        assert "reflink" not in systems["2.A.1.1.1"]

    def test_sidebar_ids_excluded(self):
        _, systems = m.parse_family_page(FAMILY_PAGE_FIXTURE)
        assert "9.Z.9.9.9" not in systems

    def test_family_filter_drops_off_family_ids(self):
        subfams, systems = m.parse_family_page(FAMILY_PAGE_FIXTURE, family="2.A.1")
        assert set(systems) == {"2.A.1.1.1"}
        assert set(subfams) == {"2.A.1.1", "2.A.1.2"}
        # filtering to a different family drops everything
        subfams9, systems9 = m.parse_family_page(FAMILY_PAGE_FIXTURE, family="9.B.1")
        assert subfams9 == {} and systems9 == {}

    def test_unnamed_family_has_no_headers_but_systems_parse(self):
        subfams, systems = m.parse_family_page(UNNAMED_FAMILY_PAGE_FIXTURE)
        assert subfams == {}
        assert systems["1.A.11.1.1"] == "Ammonia transporter and regulatory sensor, AmtB."
