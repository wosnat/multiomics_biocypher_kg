# SignalP 6.0 — Installation

SignalP 6.0 predicts signal peptides (Sec/SPI, lipoprotein, Tat) in protein
sequences. It runs inside a Docker container to avoid Python-version conflicts
(requires Python ≤3.10, PyTorch <2.0).

The `Dockerfile` in this directory mirrors the [official installation
steps](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md):
`pip install signalp-6-package/`, pin `numpy<2`, then copy the bundled model
weights into the installed `signalp` module. The built image lives in
Docker's daemon storage (shared across projects); the licensed tarball and
build context live under `~/tools/` — nothing goes into the repo.

## Prerequisites

- Docker installed and running (`docker --version` to verify)

## Step 1: Register and download the tarball

1. Go to https://services.healthtech.dtu.dk/services/SignalP-6.0/
2. Click "Downloads" (bottom of page), accept the academic license
3. Select **SignalP-6.0 fast** (the "fast" model is sufficient for bacterial proteomes)
4. DTU will email you a download link

## Step 2: Place the archive under ~/tools

```bash
mkdir -p ~/tools
mv ~/Downloads/signalp-6.0i.fast.tar.gz ~/tools/
```

(If your version differs, set `SIGNALP_TARBALL=/path/to/file.tar.gz` when
running `build.sh`.)

## Step 3: Build the Docker image

```bash
bash .claude/skills/signalp-run/build.sh
```

This:

- Stages `~/tools/signalp-6.0i/` as the build context (creates the directory
  on first run; safe to re-run — files are overwritten)
- Copies the tarball in as `signalp.tar.gz` and the skill's `Dockerfile`
- Runs `docker build -t signalp6 ~/tools/signalp-6.0i/`

First build takes 3–10 minutes (downloads Python 3.10 + PyTorch <2.0
~1 GB). Subsequent rebuilds reuse Docker's cache.

The final image (`signalp6`, ~5 GB) lives in Docker's daemon storage; the
build context (~1.5 GB) stays under `~/tools/signalp-6.0i/`.

## Step 4: Verify

```bash
docker run --rm signalp6 --help
```

You should see SignalP's help with `--fastafile`, `--organism`, `--mode`.

## Step 5: Smoke-test on a single strain

```bash
uv run python .claude/skills/signalp-run/run_signalp.py --strain HOT1A3
```

Results land in `cache/data/Alteromonas/genomes/HOT1A3/signalp/prediction_results.txt`.

## Step 6: Run on all strains

```bash
uv run python .claude/skills/signalp-run/run_signalp.py
```

## GPU upgrade path (NVIDIA hosts only)

This host has no NVIDIA GPU, so the image is CPU-only. SignalP's
[`signalp6_convert_models gpu`](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#converting-to-gpu)
step is CUDA-only and will fail without an NVIDIA driver + `nvidia-container-runtime`.

To produce a GPU image on a CUDA host:

1. Change `FROM python:3.10-slim` in `Dockerfile` to a CUDA base, e.g.
   `nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04` (then install Python 3.10).
2. Pin a CUDA-enabled PyTorch (<2.0), e.g.
   `pip install torch==1.13.1+cu117 --index-url https://download.pytorch.org/whl/cu117`.
3. After the model-weights copy step, append:
   `RUN signalp6_convert_models gpu`.
4. Tag the result separately: `SIGNALP_IMAGE_TAG=signalp6-gpu bash build.sh`.
5. Pass `--gpus all` to `docker run` (and update `run_signalp.py` accordingly).

## Troubleshooting

**"permission denied" on output files**: `run_signalp.py` already passes
`--user $(id -u):$(id -g)`. If you run `docker run` directly, add the same flag.

**"signalp tarball not found"**: `build.sh` looks for `~/tools/signalp-6.0i.fast.tar.gz`
by default. Override with `SIGNALP_TARBALL=/path/to/file.tar.gz`.

**Image size**: ~5 GB total (PyTorch dominates). Check with `docker images signalp6`.
