# Docker Cheatsheet

## Full Pipeline (build + import + post-process + deploy)

```bash
docker compose up -d
```

## Rebuild from Scratch

```bash
docker compose down -v          # -v removes the named volume (wipes the DB)
docker compose up -d
```

## Run Individual Stages

```bash
docker compose up build             # only build (create CSVs)
docker compose up import            # build + import into Neo4j
docker compose up post-process      # build + import + post-process
docker compose up deploy            # full pipeline + expose Neo4j
```

## View Logs

```bash
docker compose logs build
docker compose logs import
docker compose logs post-process    # shows homolog/expression edge counts
docker compose logs deploy
docker compose logs -f deploy       # follow live logs
```

## Access Neo4j

- Browser: http://localhost:7474
- Bolt: bolt://localhost:7687
- Auth: disabled (no username/password needed)

## Quick Cypher via CLI

```bash
docker exec deploy cypher-shell "MATCH (n) RETURN labels(n)[0] AS label, count(*) ORDER BY count(*) DESC;"
```

## Re-run Only Post-Process

If you changed `scripts/post-import.cypher` and want to re-run without rebuilding:

```bash
docker compose stop deploy
docker compose rm -f post-process
docker compose up post-process
docker compose up -d deploy
```

## Check Import Report

```bash
docker compose exec deploy cat /data/build2neo/import.report
# or locally:
cat output/import.report
```

## Inspect the Volume

```bash
docker volume inspect multiomics_biocypher_kg_biocypher_neo4j_volume
```

## Biochatter UI

- URL: http://localhost:8501
- Requires `OPENAI_API_KEY` in `.env`

## Common Issues

| Problem | Fix |
|---------|-----|
| Post-process created 0 edges | Check `docker compose logs post-process` for counts and cypher-shell errors |
| Port 7474 already in use | `docker compose down` or stop the conflicting container |
| Stale data after code changes | `docker compose down -v && docker compose up -d` to rebuild from scratch |
| Build fails (missing data) | Check that `data/` and `cache/` directories have the required files |
