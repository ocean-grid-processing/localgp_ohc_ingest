# ohc_ingest — pure-Rust, no C/system deps, so both stages stay slim.

# ---- build stage ----
FROM rust:1-slim-bookworm AS build
WORKDIR /app
COPY Cargo.toml ./
COPY src ./src
# (no Cargo.lock committed yet; cargo will resolve on first build)
RUN cargo build --release

# ---- runtime stage ----
FROM debian:bookworm-slim
COPY --from=build /app/target/release/ohc_ingest /usr/local/bin/ohc_ingest
WORKDIR /work
ENTRYPOINT ["ohc_ingest"]

# Build:
#   docker image build -t ohc_ingest ohc_pipeline/ohc_ingest
#
# Run a single month (Aug 2016) of one layer, mounting data + output:
#   docker container run --rm \
#     -v /host/FullField:/in_mean:ro \
#     -v /host/FullFieldLocalCondSim:/in_ens:ro \
#     -v /host/data:/aux:ro \
#     -v /host/out:/out \
#     -e OHC_DIR_MEAN=/in_mean -e OHC_DIR_ENSEMBLE=/in_ens -e OHC_DIR_OUT=/out \
#     -e OHC_ETOPO=/aux/etopo60.cdf -e OHC_BASINMASK=/aux/basinmask_04.msk \
#     ohc_ingest --tag OP20260110 --layer 15-20 --years 2016 --months 8
#
# Or with a config file:
#   docker run --rm -v $PWD:/work -v ... ohc_ingest config.toml --months 8
