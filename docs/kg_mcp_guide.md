# Multiomics KG — Alpha-tester MCP guide

> **Audience:** alpha testers of the *Prochlorococcus* / *Alteromonas* multi-omics
> knowledge graph who want to drive the KG from their own LLM agent (Claude Code,
> Claude Desktop, or any MCP-compatible client) via the **explorer MCP server**.
>
> **Status:** alpha. The KG content, the explorer MCP tool set, and the connection
> method may change between releases. Use the compatibility check (§5) before every
> session.

---

## 1. What this is

A live Neo4j knowledge graph integrating:

- ~80K genes across 25 cyanobacterial + heterotroph strains (+ 5 treatment / 2 reference-proteome organisms),
- ~227K differential-expression edges from ~35 publications (RNA-seq, proteomics, metabolomics, microarray),
- functional annotations (GO, KEGG, EC, COG, Pfam, TCDB, CAZy, BRITE), ortholog groups, metabolism (~3K reactions, ~1.5K metabolites), gene clusters, derived metrics, and metabolite assays.

You query it through the **explorer MCP server** (`multiomics-kg-mcp`), which exposes ~37 domain-specific tools (gene resolution, expression lookups, ontology enrichment, metabolite searches, …) on top of the raw Neo4j graph.

The graph is read-only for you. The MCP rejects every write attempt.

---

## 2. Connection details for this release

> **This is the only section that changes per host.** Everything else in this guide
> stays the same whether the KG is served from a local lab box (Track A) or from
> Neo4j Aura (Track B). The per-release values (URI, credentials, contact) come from
> the **GitHub Release notes** for the version you were invited to — start there:
> <https://github.com/wosnat/multiomics_biocypher_kg/releases>.

### 2.1 Bolt URI

| Host | URI shape | Example | TLS |
|---|---|---|---|
| **Local lab box** (Track A) | `bolt://<lab-ip>:<port>` | `bolt://132.75.249.47:17687` | no (trusted lab subnet) |
| **Aura** (Track B) | `neo4j+s://<instance-id>.databases.neo4j.io` | `neo4j+s://abc12345.databases.neo4j.io` | yes (built-in) |

Use the URI exactly as printed in the release notes; the scheme prefix (`bolt://` vs `neo4j+s://`) carries the TLS contract.

### 2.2 Reach the host (before installing anything)

Run the appropriate one-liner from the machine you'll be querying *from*:

- **Linux / macOS:** `nc -zv <host> <port>` → expect `Connection to <host> <port> port [tcp/*] succeeded!`
- **Windows 11 / 10:** `Test-NetConnection -ComputerName <host> -Port <port>` — run in **PowerShell**, not Command Prompt — expect `TcpTestSucceeded : True`.

If the test fails:

| Host | Likely cause |
|---|---|
| Local lab box | You're not on the lab subnet. Go to the lab, or connect to the university VPN (which places you back on the campus network). Off-subnet hosts are dropped at the host firewall. |
| Aura | Some campus / corporate firewalls block outbound 7687. Try from a different network, or ask your IT to whitelist `neo4j+s://*.databases.neo4j.io`. |

### 2.3 Credentials

| Host | Model | Where to get them |
|---|---|---|
| **Local lab box** | One shared `explorer` read login (not the `neo4j` admin — that's operator-only). | Distributed out-of-band by the operator (lab Slack DM, not committed to any public location). |
| **Aura** | Per-user account on the Aura instance. | Username + password from the operator; can be rotated per-tester. |

The `NEO4J_DATABASE` value is `neo4j` for both hosts unless the release notes say otherwise.

### 2.4 Read-only — what enforces it

The explorer MCP rejects every write attempt — `run_cypher` blocks any `CREATE / MERGE / SET / DELETE / DROP / CALL <write-procedure>`, and the ~36 typed tools never issue writes.

| Host | Second-layer enforcement |
|---|---|
| **Local lab box** (Neo4j Community) | **None at the DB layer** — Community has no roles, so the shared `explorer` credential is full-access if used outside the MCP. Stay in the MCP. |
| **Aura** | Your account has the Neo4j `reader` role, so writes are also rejected server-side. |

In both cases the graph is rebuilt from public sources on every release, so even a successful unauthorized write would have no lasting effect.

---

## 3. Install the MCP

The explorer MCP is published in the public repo `wosnat/multiomics_explorer`. Until we host it remotely, run it locally with `uv`:

```bash
# Requires Python 3.11+ and uv (https://docs.astral.sh/uv/).
uvx --from git+https://github.com/wosnat/multiomics_explorer multiomics-kg-mcp --help
```

If `--help` runs cleanly, you're set.

---

## 4. Configure your MCP client

The MCP reads connection settings from environment variables. Drop this into your client's MCP config (Claude Code: `.mcp.json` in the project root; Claude Desktop: `claude_desktop_config.json`). Fill the four placeholders from **§2** of this guide + the release notes.

```json
{
  "mcpServers": {
    "multiomics-kg": {
      "command": "uvx",
      "args": [
        "--from",
        "git+https://github.com/wosnat/multiomics_explorer",
        "multiomics-kg-mcp"
      ],
      "env": {
        "NEO4J_URI":      "<from §2.1 — bolt://... or neo4j+s://...>",
        "NEO4J_USERNAME": "<from §2.3>",
        "NEO4J_PASSWORD": "<from §2.3>",
        "NEO4J_DATABASE": "neo4j"
      }
    }
  }
}
```

Restart your MCP client. The `multiomics-kg` server appears with ~37 tools prefixed `mcp__multiomics-kg__…`.

---

## 5. First query — compatibility check

Before doing anything else, ask your agent to run:

```
Show me the Schema_info node — version, mcp_min_version, built_at, and the counts.
```

The agent should hit `mcp__multiomics-kg__run_cypher` with:

```cypher
MATCH (s:Schema_info {id: 'schema_info'})
RETURN s.version            AS kg_version,
       s.built_at           AS built_at,
       s.git_sha_short      AS git_sha,
       s.mcp_min_version    AS mcp_min_version,
       s.paper_count        AS papers,
       s.experiment_count   AS experiments,
       s.gene_count         AS genes,
       s.organism_count     AS organisms,
       s.expression_edge_count AS expression_edges,
       s.release_notes_url  AS release_notes
```

Confirm:
- `kg_version` matches the alpha you were invited to (NOT `0.0.0-dev` — that's a dev build, not a release).
- The installed explorer-MCP version is ≥ `mcp_min_version`. If it's older, `uvx` will pick up the latest from `main` on next launch; pin to a tag for reproducibility (§9).
- The counts roughly match the release notes — wildly different numbers mean something didn't restore correctly.

---

## 6. Find your way around

The MCP ships its own documentation as resources. Ask your agent to fetch them; they're more reliable than this doc because they ship with the MCP version you installed.

| Resource | Purpose |
|---|---|
| `docs://guide/start_here` | Decision tree: pick the right tool for your question |
| `docs://guide/concepts` | KG data model (Gene, Experiment, DerivedMetric, Metabolite, …) |
| `docs://guide/conventions` | Cross-tool semantics: `not_found` vs `not_matched`, tested-absent rows, exclude-wins-on-overlap, rankable-gated filters, AQ / informative_only defaults |
| `docs://guide/python_api` | Scripting the Python package directly (not strictly MCP, but useful for batch work) |
| `docs://tools/{tool_name}` | Per-tool reference (e.g. `docs://tools/differential_expression_by_gene`) |
| `docs://analysis/{name}` | Worked analyses (e.g. `docs://analysis/enrichment`, `docs://analysis/metabolites`) |
| `docs://examples/{file}` | Runnable example scripts |

A good opening prompt for your agent:

> Read `docs://guide/start_here` and `docs://guide/conventions`, then summarise what tools I'd use to answer "which Prochlorococcus MED4 genes change expression under nitrogen limitation, and what pathways are they in?"

---

## 7. Starter queries

Examples to verify your connection is healthy and to feel out the graph shape.

**Gene lookup:**
> Resolve `PMM0001` and show me its overview, including how many expression edges and metabolites it touches.

**Differential expression:**
> Find genes in *Prochlorococcus* MED4 whose expression changes under nitrogen limitation across all available experiments, ranked by |log2FoldChange|.

**Cross-organism orthology:**
> Take the *Prochlorococcus* MED4 katG gene and show its orthologs in other strains — does the gene exist at all in Prochlorococcus, and is it present in any of the heterotrophs in the graph?

**Metabolite search:**
> Find all metabolites measured experimentally in any *Alteromonas* strain and report which compartments they were measured in.

**Pathway enrichment:**
> I have these 50 gene locus tags — run a KEGG-pathway enrichment using the appropriate organism background.

---

## 8. Conventions and gotchas

The full conventions live at `docs://guide/conventions`. Highlights:

- **Labels are PascalCase** in Cypher: `Gene`, `Experiment`, `OrthologGroup`, `Schema_info`. Relationship types: `Changes_expression_of`, `Has_experiment`, `Gene_in_ortholog_group`, etc.
- **Booleans are strings.** `significant`, `is_time_course`, `rankable`, `git_dirty` are `"true"` / `"false"` strings, not booleans (Neo4j-array/null interaction).
- **Time-course experiments** carry parallel arrays on the `Experiment` node (`time_point_labels`, `time_point_orders`, `time_point_hours`, `time_point_totals`, …). Sentinel values: `""` = no label, `-1.0` = unknown hours.
- **`adjusted_p_value` can be null** on `Changes_expression_of` edges when the original study didn't report it. Filter accordingly.
- **`Schema_info.schema_info`** is a JSON string containing the full schema. Parse it with `apoc.convert.fromJsonMap` or just have your agent read it for richer schema introspection than `kg_schema` provides.

---

## 9. Limits and support

- **Read-only.** Enforced by the MCP for both hosts; additionally enforced server-side on Aura via the `reader` role (§2.4).
- **Query timeout.** The MCP applies a 60-second default. On Aura, the server applies its own limit too. Heavy unbounded queries get cancelled — narrow with a label, organism, or `LIMIT`.
- **Result cap.** `run_cypher` defaults to 25 rows; pass `limit=N` to bump it.
- **Bug reports / questions.** File an issue at <https://github.com/wosnat/multiomics_biocypher_kg/issues> and include the output of the compatibility check (§5) so we know which release you hit.

---

## 10. Updating

Releases are tracked on GitHub: <https://github.com/wosnat/multiomics_biocypher_kg/releases>. Each `kg-X.Y.Z[-(alpha|beta|rc).N]` tag has a release page with the changelog entry, build metadata (`metadata.json`), and the **per-release connection values** described in §2. **Watch the repo → Custom → Releases** to get email notifications.

You don't download anything — the data lives on the host (§2). When a new release ships, the release page says whether the URI / credentials have changed; update your `.mcp.json` accordingly.

When you re-run `uvx --from git+...` it pulls the latest explorer-MCP commit; pin to a tag with `@kg-0.1.0-alpha.1` if you need reproducibility.

Re-run the compatibility check (§5) after every release notification — schema may have moved.
