# Build the guinea_nofftw binary in a throwaway build stage.
FROM ubuntu:22.04 AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc make libc6-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY . .
RUN make guinea_nofftw

# Minimal runtime image with everything pre-installed and configured so a
# plain `docker run` produces simulated events with no further setup.
FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends python3 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /build/guinea_nofftw ./guinea_nofftw
COPY --from=builder /build/acc.dat ./acc.dat
COPY --from=builder /build/test_params.dat ./test_params.dat
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
COPY docker/make_hepmc.py /usr/local/bin/make_hepmc.py
RUN chmod +x /usr/local/bin/entrypoint.sh /usr/local/bin/make_hepmc.py

ENV ACCELERATOR=muon \
    PARAMS=muon_pairs_10tev \
    N_EVENTS=1 \
    OUTPUT_DIR=/output

VOLUME ["/output"]

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
