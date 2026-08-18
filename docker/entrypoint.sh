#!/bin/bash
# Entrypoint for the guineapig_mumu Docker image.
#
# Runs the GuineaPig 10 TeV muon pair-production simulation and writes the
# resulting log and output file(s) both to stdout and to a mountable output
# directory, with no further action required from the user.
#
# In addition to the raw GuineaPig output, a HepMC3 ASCII event file is
# produced for every bunch crossing by converting the beam energies and the
# stored photons/pair-production leptons into a HepMC event record. By
# default, the per-event HepMC files for the run are also merged into a
# single multi-event HepMC file.
#
# Configuration (env vars, all optional):
#   ACCELERATOR  - accelerator definition from acc.dat (default: muon)
#   PARAMS       - parameter set from acc.dat/test_params.dat (default: muon_pairs_10tev)
#   N_EVENTS     - number of bunch crossings ("events") to simulate (default: 1)
#   OUTPUT_DIR   - directory the log and output files are written to (default: /output)
#   PT_MIN       - pt cut in GeV applied to the charged pair leptons during the
#                  HepMC conversion (default: 0.015; MAIA inside-beam-pipe
#                  value is 0.017). Set to 0 to disable.
#   MERGE_HEPMC  - whether to merge the per-event HepMC files for the run into
#                  a single "<params>.hepmc" file (default: 1). Set to 0 to
#                  skip merging and keep only the per-event files.
#
# N_EVENTS may also be given as the first positional argument, e.g.:
#   docker run ghcr.io/<owner>/guineapig_mumu 5
set -euo pipefail

ACCELERATOR="${ACCELERATOR:-muon}"
PARAMS="${PARAMS:-muon_pairs_10tev}"
N_EVENTS="${N_EVENTS:-1}"
OUTPUT_DIR="${OUTPUT_DIR:-/output}"
PT_MIN="${PT_MIN:-0.015}"
MERGE_HEPMC="${MERGE_HEPMC:-1}"

if [ $# -ge 1 ]; then
    N_EVENTS="$1"
fi

if ! [[ "$N_EVENTS" =~ ^[0-9]+$ ]] || [ "$N_EVENTS" -lt 1 ]; then
    echo "N_EVENTS must be a positive integer (got: $N_EVENTS)" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

WORKDIR="$(mktemp -d)"
cp /app/guinea_nofftw /app/acc.dat /app/test_params.dat "$WORKDIR/"
cat "$WORKDIR/test_params.dat" >> "$WORKDIR/acc.dat"
cd "$WORKDIR"

LOG_FILE="$OUTPUT_DIR/guinea_pig.log"
: > "$LOG_FILE"

echo "Simulating $N_EVENTS bunch crossing(s) with accelerator='$ACCELERATOR' parameters='$PARAMS'" | tee -a "$LOG_FILE"

for i in $(seq 1 "$N_EVENTS"); do
    OUT_NAME="${PARAMS}_event${i}.out"

    {
        echo
        echo "=== Event $i/$N_EVENTS ==="
    } | tee -a "$LOG_FILE"

    ./guinea_nofftw "$ACCELERATOR" "$PARAMS" "$OUT_NAME" 2>&1 | tee -a "$LOG_FILE"

    cp "$OUT_NAME" "$OUTPUT_DIR/$OUT_NAME"

    for extra in photon.dat lumi.ee.out lumi.eg.out lumi.ge.out lumi.gg.out pairs.dat pairs0.dat; do
        if [ -f "$extra" ]; then
            cp "$extra" "$OUTPUT_DIR/event${i}_${extra}"
        fi
    done

    HEPMC_NAME="${PARAMS}_event${i}.hepmc"
    PAIRS_FILE=""
    if [ -f pairs.dat ]; then
        PAIRS_FILE="pairs.dat"
    elif [ -f pairs0.dat ]; then
        PAIRS_FILE="pairs0.dat"
    fi

    python3 /usr/local/bin/make_hepmc.py \
        --out "$OUT_NAME" \
        --photons photon.dat \
        ${PAIRS_FILE:+--pairs "$PAIRS_FILE"} \
        --event-number "$i" \
        --pt-min "$PT_MIN" \
        --output "$OUTPUT_DIR/$HEPMC_NAME" 2>&1 | tee -a "$LOG_FILE"

    {
        echo
        echo "--- Output file: $OUT_NAME ---"
        cat "$OUT_NAME"
    } | tee -a "$LOG_FILE"
done

if [[ "$MERGE_HEPMC" != "0" && "$N_EVENTS" -gt 0 ]]; then
    MERGED_HEPMC="$OUTPUT_DIR/${PARAMS}.hepmc"
    python3 /usr/local/bin/merge_hepmc.py \
        --output "$MERGED_HEPMC" \
        "$OUTPUT_DIR/${PARAMS}"_event*.hepmc 2>&1 | tee -a "$LOG_FILE"
fi

echo | tee -a "$LOG_FILE"
echo "Done. Log and output file(s) written to $OUTPUT_DIR" | tee -a "$LOG_FILE"
