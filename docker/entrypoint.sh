#!/bin/bash
# Entrypoint for the guineapig_mumu Docker image.
#
# Runs the GuineaPig 10 TeV muon pair-production simulation and writes the
# resulting log and output file(s) both to stdout and to a mountable output
# directory, with no further action required from the user.
#
# In addition to the raw GuineaPig output, a HepMC2 ASCII event file is
# produced for every bunch crossing by converting the beam energies and the
# stored photons/pair-production leptons into a HepMC event record.
#
# Configuration (env vars, all optional):
#   ACCELERATOR  - accelerator definition from acc.dat (default: muon)
#   PARAMS       - parameter set from acc.dat/test_params.dat (default: muon_pairs_10tev)
#   N_EVENTS     - number of bunch crossings ("events") to simulate (default: 1)
#   OUTPUT_DIR   - directory the log and output files are written to (default: /output)
#   PT_MIN       - pt cut in GeV applied to the charged pair leptons during the
#                  HepMC conversion (default: 0.015; MAIA inside-beam-pipe
#                  value is 0.017). Set to 0 to disable.
#
# N_EVENTS may also be given as the first positional argument, e.g.:
#   docker run ghcr.io/<owner>/guineapig_mumu 5
set -euo pipefail

ACCELERATOR="${ACCELERATOR:-muon}"
PARAMS="${PARAMS:-muon_pairs_10tev}"
N_EVENTS="${N_EVENTS:-1}"
OUTPUT_DIR="${OUTPUT_DIR:-/output}"
PT_MIN="${PT_MIN:-0.015}"

if [ $# -ge 1 ]; then
    N_EVENTS="$1"
fi

if ! [[ "$N_EVENTS" =~ ^[0-9]+$ ]] || [ "$N_EVENTS" -lt 1 ]; then
    echo "N_EVENTS must be a positive integer (got: $N_EVENTS)" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

WORKDIR="$(mktemp -d)"
cp /app/guinea /app/acc.dat /app/test_params.dat "$WORKDIR/"
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

    ./guinea "$ACCELERATOR" "$PARAMS" "$OUT_NAME" 2>&1 | tee -a "$LOG_FILE"

    cp "$OUT_NAME" "$OUTPUT_DIR/$OUT_NAME"

    for extra in photon.dat lumi.ee.out lumi.eg.out lumi.ge.out lumi.gg.out pairs.dat pairs0.dat; do
        if [ -f "$extra" ]; then
            cp "$extra" "$OUTPUT_DIR/event${i}_${extra}"
        fi
    done

    HEPMC_NAME="${PARAMS}_event${i}.hepmc"
    # pairs.dat is the TRACKED output (print_pairs(), needs track_pairs>0).
    # pairs0.dat is written at production time, before the opposing beam's field
    # deflects the pairs, and that deflection is what sets their pT -- using it
    # understates detector occupancy by orders of magnitude. Falling back to it
    # silently is how that went unnoticed, so the fallback is now loud.
    PAIRS_FILE=""
    if [ -f pairs.dat ]; then
        PAIRS_FILE="pairs.dat"
    elif [ -f pairs0.dat ]; then
        PAIRS_FILE="pairs0.dat"
        {
            echo
            echo "############################################################"
            echo "WARNING: pairs.dat not found -- falling back to pairs0.dat."
            echo "  pairs0.dat holds PRODUCTION-TIME kinematics with no beam-field"
            echo "  deflection (mean pT 1.7 MeV vs 39 MeV tracked). The resulting"
            echo "  HepMC is NOT suitable for detector background studies."
            echo "  Set track_pairs=1 (and pair_step=5.0) in the parameter set."
            echo "############################################################"
            echo
        } | tee -a "$LOG_FILE" >&2
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

echo | tee -a "$LOG_FILE"
echo "Done. Log and output file(s) written to $OUTPUT_DIR" | tee -a "$LOG_FILE"
