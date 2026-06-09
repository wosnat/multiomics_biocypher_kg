# Knowledge Graph Alpha — IT Review & Approval Request

**Prepared by:** «your name», «lab / department»
**Contact:** «email»
**Date:** 2026-05-24
**Machine:** «hostname» — `132.75.249.47` (existing lab workstation, already on the campus network)

---

## 1. Summary

We would like to run a **read-only graph database** (Neo4j, in Docker) on our existing lab
workstation so that **~4–5 members of our group** can query a scientific knowledge graph from
their own analysis tools during an internal alpha test.

- It is **internal and small** (≤5 users, all in our group).
- It contains **only public scientific data** — no personal, human-subjects, confidential, or regulated data.
- It is **not exposed to the public internet**. Access is restricted to our **lab subnet** only.
- It requires **no changes at the campus perimeter firewall** and **no new inbound rules from outside the university**.

We are asking for confirmation that this is acceptable under university IT/security policy, plus a
couple of small networking conveniences (a stable IP and the lab subnet range). Details below.

---

## 2. What we want to run

| | |
|---|---|
| **Software** | Neo4j 5.15 **Community Edition**, official Docker image, plus our own data-loading pipeline |
| **Where** | The existing lab workstation `132.75.249.47` (no new hardware) |
| **What's new** | A second, **read-only** Neo4j instance ("alpha") alongside the existing development instance, which stays unchanged and **localhost-only** |
| **Purpose** | Let group members query a multi-omics knowledge graph during an internal alpha |
| **Duration** | Ongoing internal alpha; can be stopped at any time |

The existing development instance is bound to `localhost` only and is **not** part of this request — nothing about it changes.

---

## 3. Data sensitivity — none

The graph is built entirely from **publicly available scientific data**: genome, gene-expression,
proteomics, and metabolomics data for marine bacteria (*Prochlorococcus*, *Alteromonas*),
assembled from public reference databases (NCBI, UniProt, KEGG and similar) and from published,
peer-reviewed papers.

There is **no personal data, no human-subjects data, no confidential or export-controlled data,
and no credentials or secrets** in the graph. The entire dataset is already public, and the graph
is **rebuilt from those public sources** for each release — nothing of record exists only here.

---

## 4. Network exposure

| Item | Detail |
|---|---|
| **Ports** | `17687` (Neo4j Bolt, the query protocol), `17474` (Neo4j HTTP browser) — deliberately **non-default** ports |
| **Bind address** | The machine's **lab IP only** (`132.75.249.47`) — **not** `0.0.0.0`, not all interfaces |
| **Reachable from** | The **lab subnet only** (`«lab subnet CIDR, e.g. 132.75.249.0/24»`), enforced by a host firewall rule on the machine |
| **Public internet** | **No exposure.** No inbound perimeter rule is requested; the campus firewall already blocks external access |
| **Future (optional)** | If group members later need access from home, they would use the **existing university VPN** (which places them back on the campus network). This would be enabled simply by adding the VPN address range to the same host firewall allowlist — **still no perimeter change** |

In short: the service answers only to hosts already inside our lab network. It is not, and will not
be, reachable from the public internet.

---

## 5. Access control

- **Authentication is enabled** — Neo4j username/password. No anonymous access.
- **Read-only for users.** Group members connect with a **shared read credential** (distributed privately, not stored in any public location). Day-to-day access is read-only, enforced by the query tool they use.
- **Low write-risk by design.** Because the graph holds only public data and is **rebuilt from source** on each release, any accidental or malicious change is **low-impact and fully recoverable** — there is no unique or sensitive data to protect.
- The **administrative credential is held only by the operator** («your name») and is not shared with users.

---

## 6. Security posture (why this is low-risk)

- **Small, internal, read-only**: ≤5 known group members; public data; no writes expected.
- **Not internet-facing**: bound to the lab IP and firewalled to the lab subnet; no perimeter changes.
- **Official, current image**: Neo4j 5.15 Community from the official Docker registry; we can apply updates. Attack surface is minimal because the service is not internet-exposed and holds no sensitive data.
- **Non-default ports** (`17687`/`17474`) and **authentication on**.
- **Recoverable**: the graph is regenerated from public sources; loss or tampering has no lasting impact.
- **Clear ownership**: maintained by «your name», «email».

---

## 7. What we are asking IT to approve / provide

1. **Approval** to run this read-only Neo4j (Docker) service on `132.75.249.47`, reachable on the **lab subnet** on ports `17687` and `17474`.
2. A **static IP or DHCP reservation** for «hostname» so its address is stable for users.
3. **Confirmation of the lab subnet range (CIDR)** so we scope the host firewall allowlist correctly (we want to restrict to the lab, not the whole campus).
4. **Confirmation** that the above is consistent with university IT / security / acceptable-use policy.

We're happy to adjust any of this to meet your requirements — e.g. tighten the allowed source
range, change ports, or add TLS — and to walk through it in person.

---

*Internal technical detail (build pipeline, release process, schema) is documented separately in
`plans/alpha_release.md`; this document covers only what is relevant to IT review.*
