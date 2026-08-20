# Build the guinea_nofftw binary in a throwaway build stage.
FROM ubuntu:22.04 AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc make libc6-dev libfftw3-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY . .
# Build the FFTW target: guinea_nofftw omits fourtrans3.c (the FFT Poisson
# solver), so the two builds compute different fields and must not be mixed.
RUN make guinea FFTW_HOME=/usr

# Minimal runtime image with everything pre-installed and configured so a
# plain `docker run` produces simulated events with no further setup.
FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends python3 libfftw3-3 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /build/guinea ./guinea
COPY --from=builder /build/acc.dat ./acc.dat
COPY --from=builder /build/test_params.dat ./test_params.dat
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
COPY docker/make_hepmc.py /usr/local/bin/make_hepmc.py
RUN chmod +x /usr/local/bin/entrypoint.sh /usr/local/bin/make_hepmc.py

ENV ACCELERATOR=muon \
    PARAMS=muon_pairs_10tev \
    N_EVENTS=1 \
    PT_MIN=0.015 \
    OUTPUT_DIR=/output

VOLUME ["/output"]

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
