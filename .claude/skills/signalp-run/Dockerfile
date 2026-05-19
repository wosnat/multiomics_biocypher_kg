FROM python:3.10-slim

# SignalP 6.0 (fast) — installed per official instructions:
# https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md
#
# build.sh renames the versioned upstream archive (e.g.
# signalp-6.0i.fast.tar.gz) to `signalp.tar.gz` inside the build context so
# this Dockerfile stays version-agnostic.
COPY signalp.tar.gz /tmp/

RUN cd /tmp && tar xzf signalp.tar.gz \
    && PKG_DIR=$(find /tmp -maxdepth 3 -type d -name signalp-6-package | head -1) \
    && test -n "$PKG_DIR" || (echo "ERROR: signalp-6-package not found in tarball" && exit 1) \
    && pip install --no-cache-dir "numpy<2" \
    && pip install --no-cache-dir "$PKG_DIR/" \
    && SIGNALP_DIR=$(python -c "import signalp, os; print(os.path.dirname(signalp.__file__))") \
    && cp -r "$PKG_DIR/models/"* "$SIGNALP_DIR/model_weights/" \
    && rm -rf /tmp/signalp* /tmp/signalp.tar.gz

WORKDIR /data

ENTRYPOINT ["signalp6"]
