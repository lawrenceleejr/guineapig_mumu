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
#   PT_MIN       - extra flat pt cut in GeV on the charged pair leptons during
#                  the HepMC conversion (default: 0, i.e. off -- the geometric
#                  filter below supersedes it).
#   GEOM_FILTER  - drop pair leptons whose helix in the solenoid field can
#                  never reach the bore wall inside the detector, and so can
#                  produce neither a hit nor a shower (default: 1). This is
#                  loss-free for the bore described by BORE_*; set to 0 to
#                  keep every lepton.
#   B_FIELD      - solenoid field in T for that filter (default: 5, MAIA).
#   BORE_R0      - bore inner radius at the IP in mm (default: 23).
#   BORE_SLOPE   - rate the bore opens with |z|, mm/m (default: 10).
#   BORE_ZMAX    - |z| in m at which a particle has left the region (default: 6).
#                  BORE_* MUST match the geometry being simulated; a real bore
#                  narrower than BORE_R0 would make the filter lossy.
#   DROP_PHOTONS - omit beamstrahlung photons, which are 0.3% of the record and
#                  all within 3 mrad of the axis (default: 0, i.e. keep them).
#   HEPMC_WEIGHT - HepMC event weight. Left unset it is derived from the
#                  pair_ratio GuineaPig reports, so a thinned run carries its
#                  own 1/pair_ratio scale factor.
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
PT_MIN="${PT_MIN:-0}"
MERGE_HEPMC="${MERGE_HEPMC:-1}"
GEOM_FILTER="${GEOM_FILTER:-1}"
B_FIELD="${B_FIELD:-5.0}"
BORE_R0="${BORE_R0:-23.0}"
BORE_SLOPE="${BORE_SLOPE:-10.0}"
BORE_ZMAX="${BORE_ZMAX:-6.0}"
DROP_PHOTONS="${DROP_PHOTONS:-0}"

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

    # A pair_ratio < 1 run stores only that fraction of the pairs, so each
    # stored lepton stands for 1/pair_ratio real ones. Read the value back out
    # of GuineaPig's own output rather than trusting the caller to remember.
    if [ -z "${HEPMC_WEIGHT:-}" ]; then
        PAIR_RATIO="$(sed -n 's/.*pair_ratio=\([0-9.eE+-]*\);.*/\1/p' "$OUT_NAME" | head -1)"
        HEPMC_WEIGHT="$(python3 -c "
import sys
r = sys.argv[1]
try:
    r = float(r)
except ValueError:
    r = 1.0
print(1.0/r if r > 0 else 1.0)" "${PAIR_RATIO:-1}")"
        if [ "$HEPMC_WEIGHT" != "1.0" ]; then
            echo "pair_ratio=$PAIR_RATIO -> HepMC event weight $HEPMC_WEIGHT" \
                | tee -a "$LOG_FILE"
        fi
    fi

    # built as an array: under `set -e` a $(cond && echo flag) that evaluates
    # false makes the substitution exit non-zero, which is a trap not worth
    # walking into for two optional flags.
    HEPMC_ARGS=(--out "$OUT_NAME" --photons photon.dat)
    [ -n "$PAIRS_FILE" ] && HEPMC_ARGS+=(--pairs "$PAIRS_FILE")
    HEPMC_ARGS+=(--event-number "$i" --pt-min "$PT_MIN"
                 --b-field "$B_FIELD" --bore-r0 "$BORE_R0"
                 --bore-slope "$BORE_SLOPE" --bore-zmax "$BORE_ZMAX"
                 --weight "$HEPMC_WEIGHT")
    [ "$GEOM_FILTER" = "0" ] && HEPMC_ARGS+=(--no-geometric-filter)
    [ "$DROP_PHOTONS" = "1" ] && HEPMC_ARGS+=(--drop-photons)
    HEPMC_ARGS+=(--output "$OUTPUT_DIR/$HEPMC_NAME")

    python3 /usr/local/bin/make_hepmc.py "${HEPMC_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE"

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
