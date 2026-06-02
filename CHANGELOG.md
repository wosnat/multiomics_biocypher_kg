# Changelog

All notable changes to the multi-omics knowledge graph are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions use the KG release scheme `X.Y.Z[-(alpha|beta|rc).N]` and are tagged
`kg-X.Y.Z…`.

**Process (accumulate-then-cut):** log notable changes under **[Unreleased]** as
they land. At release time, `/release-kg` *cuts* `[Unreleased]` into a dated
version section, stamps the same version onto `Schema_info.version`, and renders
the GitHub Release notes from that section. The changelog is the source of truth;
the GitHub Release is a rendering of one section. See `plans/alpha_release.md` §2.3.

## [Unreleased]

### Added

### Changed

### Fixed

## [0.1.0-alpha.1] - 2026-06-02

### Added
- `Schema_info` release metadata, stamped at post-import (every build): `version`,
  `built_at`, `git_sha`, `git_sha_short`, `git_branch`, `git_dirty`,
  `mcp_min_version`, `release_notes_url`, plus computed counts (`paper_count`,
  `experiment_count`, `gene_count`, `organism_count`, `expression_edge_count`).
  Dev builds stamp `0.0.0-dev`; releases stamp the tagged version. Added in both
  `scripts/post-import.sh` (Group 4) and `scripts/post-import.cypher`.

### Changed

### Fixed
