# CLAUDE.md — working memory for guineapig_mumu

Cross-session context for Claude. This repo is the **generator** end of a
10 TeV muon-collider incoherent-pair background chain; the consumer is the
**MAIA** detector simulation, which lives elsewhere (see *Cluster* below).

Read the *Verified* and *Unverified* sections before acting on any number.
Everything in *Verified* was measured in this repo; everything in
*Unverified* is an assumption that a real geometry or a real `ddsim` run can
overturn.

---

## 1. The standing goal

Shrink what GuineaPig hands to MAIA's `ddsim` **without removing anything that
would produce a detector hit**, because the simulation stage is the bottleneck.

Two separate objectives, often conflated — keep them apart:

- **Loss-free reduction.** Drop only particles that provably interact with
  nothing. Always acceptable.
- **Statistical thinning.** Drop a random fraction and scale up. Unbiased in
  the mean, wrong in the fluctuations. Acceptable for occupancy maps, *not*
  for cluster merging, tracking confusion, or occupancy tails.

---

## 2. Verified — measured in this repo, trust these

All from one tracked bunch crossing (BX) of `muon_pairs_10tev`, run with the
**FFTW `guinea` binary** (`make guinea FFTW_HOME=/usr`, links `libfftw3.so.3`).
`guinea_nofftw` computes the beam-beam field differently and must never be
mixed into a comparison.

### Composition of one BX

| | count | note |
|---|---|---|
| pair leptons (`pairs.dat`) | 560,566 | 99.7% of the record, all of the cost |
| beamstrahlung photons (`photon.dat`) | 1,810 | see below — unusable |

Pair leptons: median E **68 MeV**, median pT **9.8 MeV**, median polar angle
**6.6°**. Total pair energy **639.5 TeV**. GuineaPig runtime ≈ **59 s** per BX,
single core, so the generator is *not* the bottleneck — `ddsim` is.

The photons are useless twice over: every one measured is within **3 mrad** of
the axis (max 2.99 mrad, total 3.2 GeV) so none can reach a detector element,
**and** they are not normalised like the pairs — `store_full_pair()` samples
against the macroparticle weight `wgt`, `store_photon()` carries no weight at
all. Never quote them as a photon-background estimate.

### Cut ladder

| configuration | leptons in HepMC | file |
|---|---|---|
| no filter | 560,566 | 44.0 MB |
| flat `pT > 15 MeV` (the old default) | 238,038 | 18.1 MB |
| geometric filter, 23 mm + 10 mm/m bore | 202,537 | 15.5 MB |
| `muon_pairs_10tev_reduced` + geometric | 19,549 | 1.5 MB |

- Flat `pT > 15 MeV` keeps 42.5% of leptons but **85.6% of the energy**;
  the geometric filter keeps 36% of leptons and **38%** of the energy. The
  high-energy pairs are the collinear ones that go down the bore into the MDI.
- Filter variants on the same BX: **202,537** as shipped (conservative
  off-axis), **196,048** exact off-axis, **195,552** if the production offset
  is ignored — the last would silently lose 496 real wall-reaching particles.
- `pair_ecut=0.017` (up from 5 MeV): 415,515 stored, ×1.35, and the
  filter-surviving population moves only **0.7%** (196,048 → 194,743). Not
  exactly free: the cut is at production, before the beam-beam field does work.
- `pair_ratio=0.1`: 56,483 stored; scaled by 10 it reproduces the full sample
  to **0.1% in count, 1.4% in energy**. Exactly uniform thinning.
- `pairs.dat` columns 5–7 are position **in nm** (both transverse *and*
  longitudinal — `guinea_pig.c` prints `cut_z*1e-3` on output, so `cut_z` is
  held internally in nm too). Median 0.18 mm from the axis, max 5.6 mm.
- `pairs0.dat` and `pairs.dat` are **not row-aligned** (charge-sign agreement
  is 50% in both orderings) despite holding the same particles. Do not try to
  pair them up per-particle.

### Why the filter matters far more to ddsim than the particle count implies

| | leptons | helix turns Geant4 must propagate |
|---|---|---|
| dropped by the filter | 358,029 | **1.09 × 10⁸** |
| kept | 202,537 | 2.5 × 10⁴ |

A dropped lepton spirals a median of **43.8 turns** (p90 278) at a 2.2 mm
helix radius through 6 m of 5 T field, in vacuum, hitting nothing. The turn
count is an upper bound — Geant4's looper-killer clips the tail — but the
filter's value is orders of magnitude of field stepping, not the 2.8× in
particle count.

### Where the survivors deposit — this targets the ddsim settings

| first wall contact | leptons | energy |
|---|---|---|
| \|z\| < 6 cm | 50.5% | 2.8% |
| \|z\| < 50 cm | 92.7% | 40.5% |
| \|z\| < 2 m | 99.1% | 84.4% |

Median first contact at z = 58.8 mm after only **0.095 turns** — survivors go
nearly straight into the nozzle tip. The cost is **EM showers in tungsten in
the first half-metre**, not tracking in the detector volume.

Survivor energy: p10 23.7 MeV, median 161 MeV, p90 2.88 GeV. **40.8% are below
100 MeV but carry 1.5% of the energy** — so a kinetic-energy cut on *primaries*
is a bad trade (loses 41% of hit-producing particles, saves almost no shower
work). Cut secondaries in the nozzle instead.

### MAIA parameters used (from arXiv:2502.00181)

5 T solenoid. Vertex barrel R = 3.0–10.4 cm, |Z| ≤ 65 cm. Inner tracker
R = 12.7–55.4 cm. Outer tracker R = 81.9–148.6 cm. Detector half-length
≈ 231 cm. Nozzle shadows |η| > 2.44 (≈10° polar). Quoted BIB occupancy
≈ 30,000 hits/cm² in the innermost pixel layer.

---

## 3. Unverified — assumptions that need the real geometry

**These are the things to fix first once there is cluster access.**

1. **`BORE_R0 = 23 mm`, `BORE_SLOPE = 10 mm/m`, `BORE_ZMAX = 6 m` are not from
   the MAIA XML.** 23 mm was reverse-engineered from this repo's own 17 MeV
   figure (`pT = 0.3·B·r/2` ⇒ r = 22.7 mm). The whole filter rests on it. A
   real bore **narrower** than `BORE_R0` makes the filter lossy. Sensitivity:
   10 mm → 246,276 kept; 15 mm → 223,760; 23 mm → 196,048; 30 mm → 175,998.
2. **The bore is modelled as a straight cone** `r0 + slope·|z|`. If the real
   nozzle bore has a local minimum or a step, particles the filter clears
   could clip it.
3. **Uniform 5 T out to |z| = 6 m** is assumed. The real solenoid is bounded;
   beyond it trajectories are straight.
4. **Whether `ddsim` reads our HepMC3 at all.** `make_hepmc.py` writes HepMC3
   Asciiv3 with a `.hepmc` extension, plus a `W` event-weight line I added. If
   DD4hep's reader dispatches `.hepmc` to a HepMC2 parser it will fail or
   misparse. **Verify before trusting any downstream number.** Try
   `.hepmc3` / `.hepmc` and check the particle count round-trips.
5. **The `--doOverlayIP` input format.** The framework has an incoherent-pair
   overlay path (see below) but the docs I read only showed the file-path
   parameters for `--doOverlayFull`. The right integration is probably to
   simulate our pairs once into an EDM4hep BIB-style file and feed that, *not*
   to run our HepMC as primaries per physics event.

---

## 4. Decisions already made, and why

- **The geometric filter replaced the flat `PT_MIN` as the default.** `PT_MIN`
  still exists (default 0 = off) for reproducing old samples. The flat cut is
  the saturated limit of the same containment argument and is strictly worse.
- **The filter is closed-form, not a scan.** `f(phi)` is concave up to a half
  turn and falling after, so the maximum excursion is at one stationary point.
  One evaluation per particle; 560k in 2 s with **no numpy** — required,
  because the runtime image has none.
- **The off-axis production offset is handled conservatively** (chord must span
  `r_bore − r_prod`). Costs 3.3% extra kept particles; never drops a keeper.
  The exact treatment is not concave so it would need a per-particle scan.
- **`pair_ratio` thinning is opt-in**, via `muon_pairs_10tev_reduced`, never a
  default, because it is the one lever that is not safe by default.
- **The `1/pair_ratio` weight is written as the HepMC `W` record** and the
  entrypoint derives it from GuineaPig's own reported `pair_ratio`, so it
  cannot get lost. **Nothing applies it for you and `ddsim` ignores it** —
  hit densities must be scaled at analysis time.
- **Two upstream bug fixes must not regress** (CI guards them): `ph[10000]`
  in `step_pair_1()` (was `ph[1000]`, stack smashing), and `pair_step=5.0`
  (larger value subdivides *more* finely; 0.2 aborts with "too many photons").
  `pair_step=5.0` **has never been checked for convergence** — vary it before
  trusting absolute numbers.

### Known approximations still in the converter

- Every particle is written at the **origin**; the tracked production vertex
  (≤5.6 mm) is discarded. 0.25% effect on the filter, but it smears the
  background's true d0/z0 in the vertex detector.
- `store_pairs=2` in `muon_pairs_10tev` writes a 560k-line `pairs0.dat` that
  nothing consumes. Kept only because CI's mean-pT regression test compares
  against it. `muon_pairs_10tev_reduced` uses `store_pairs=1`.

---

## 5. Cluster / MAIA software stack

From <https://mcd-wiki.web.cern.ch/software/tutorials/fermilab2026/setup/>
(Fermilab 2026 tutorial). **v3.1 removed Marlin/ILCSoft/LCIO** — the chain is
Gaudi + EDM4hep, and `.slcio` from older productions cannot be read.

```bash
# login host (OSG access point), work dir
ssh ap23.uc.osg-htc.org
cd /scratch/$USER/tutorial2026

# container, from cvmfs on OSG
apptainer run -B /ospool/uc-shared/project/futurecolliders/data/ \
  /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24\:v3.1-amd64
# stack is auto-sourced by the entrypoint; manually: source /opt/setup_mucoll.sh
command -v k4run ddsim          # verify

# geometry
git clone --recurse-submodules --branch tutorial_20260817 \
  https://github.com/MuonColliderSoft/mucoll-benchmarks.git
cd mucoll-benchmarks && source setup_config.sh . MAIA_v0
# geometries: MAIA_v0 | MuSIC_v2 | MuColl_v1
# sets MUCOLL_GEO (compact XML), MUCOLL_GEOM_NAME, MUCOLL_CONFIG, MUCOLL_MATMAP
```

Simulation and digitisation:

```bash
ddsim --steeringFile simulation/steer_baseline.py \
      --inputFiles gen_output.edm4hep.root \
      --outputFile sim_output.edm4hep.root \
      --numberOfEvents 100
# steer_baseline.py sets the physics list, the REGION CUTS, and output collections

k4run "$MUCOLL_CONFIG/$MUCOLL_CONFIG_NAME/digi_steer.py" \
      --doOverlayFull --OverlayFullNumberBackground 10 \
      --OverlayFullPathToMuPlus  /path/to/MuPlus/files \
      --OverlayFullPathToMuMinus /path/to/MuMinus/files
# --doOverlayIP  enables INCOHERENT PAIR overlay as a separate step  <-- our sample
# --keepEverything preserves hit collections (dropped by default with overlay on)
# digitiser collection precedence: incoherent pairs > BIB > raw

k4run .../reco_steer.py --numThreads 10 --TrackingThreads 5   # needs numThreads*TrackingThreads CPUs
```

Existing 10 TeV BIB samples for `MAIA_v0`, EDM4hep:
`/ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC_MAIA_v0/v9/BIB10TeV`

Other tutorial pages: `sim/`, `digireco/`, `analysis/`, `hackathon/`. The
hackathon page notes full-BIB runs require Condor batch submission, and that
Yoke/muon collections are **not** currently overlaid (`Overlay/overlay_BIB.py`).

---

## 6. First-session checklist when cluster access arrives

In order — the first three items can invalidate work built on the rest.

1. **Extract the real bore profile.** Parse `$MUCOLL_GEO` (follow its
   `<include>`s) for the beryllium beam pipe and nozzle solids; get the inner
   radius as a function of z. Set `BORE_R0`/`BORE_SLOPE`/`BORE_ZMAX` from it,
   or replace the linear cone with a piecewise profile if the real bore is not
   conical. Re-run the sensitivity table and update this file.
2. **Confirm `ddsim` reads our HepMC.** Round-trip a small file, check the
   particle count and that the `W` line does not break the parser. If HepMC3
   is not accepted, add an EDM4hep writer to `make_hepmc.py` rather than
   downgrading to HepMC2.
3. **Inspect an existing BIB sample** under `.../DataMuC_MAIA_v0/v9/BIB10TeV`
   to learn the collection names, units and per-file BX convention that
   `--doOverlayIP` expects — then match it instead of inventing a format.
4. **Measure the actual cost**, don't trust any ranking including the one in
   §2. Time a fixed subsample (say 1,000 primaries) per setting:
   filtered vs unfiltered; then the levers below.
5. **Validate the filter empirically once.** Simulate the *dropped* population
   and confirm it produces ≈0 hits in MAIA. That is the one measurement that
   turns the filter from an argument into a fact.

### ddsim levers, ranked by the §2 measurements (all untested)

Option spellings are standard DD4hep/ddsim practice — **verify against the
container's version before wiring them in.**

1. **A `<region>` around the nozzle with a coarse range cut / secondary
   threshold**, in the compact XML or via `steer_baseline.py`'s region cuts.
   93% of primaries land within 50 cm, so this is where shower suppression
   buys almost everything and touches the tracker not at all. Preferred over
   a global `--physics.rangecut`.
2. **`G4UserLimits` on the nozzle volumes** (a `<limitset>` with a minimum
   kinetic energy) — the surgical form of the same idea.
3. **Field-propagation tolerances** (`SIM.field.delta_chord`,
   `delta_one_step`, `eps_min/max`, `largest_step`) and Geant4's
   `/process/transportation/thresholdImportantEnergy`. Normally the top BIB
   lever, but **only if not filtering** — the filter already removes the
   1.09 × 10⁸ spiral turns.
4. **Cheaper EM physics constructor.** Real speedup, degrades low-energy
   accuracy; usually fine for occupancy.
5. **Split and parallelise.** BIB is embarrassingly parallel and the merged
   multi-event HepMC splits trivially; `--numThreads`/`--TrackingThreads`
   exist on the reco side.

Do **not** reach for `SIM.part.minimalKineticEnergy` for speed — it governs
what reaches the MCParticle collection, *after* the tracking is paid for. It is
the lever for output file size (that collection is >99% of the bytes), not CPU.

---

## 7. Repo conventions

- Develop on `claude/maia-detector-output-optimization-t0s5qi`. Never push
  elsewhere without asking.
- Build with FFTW: `make guinea FFTW_HOME=/usr` (needs `libfftw3-dev`).
  Never mix `guinea` and `guinea_nofftw` results.
- `make_hepmc.py` must keep running with **stdlib only** — the runtime image
  has no numpy. Analysis scratch work may use numpy.
- CI asserts: tracked pT ≫ produced pT; the closed-form filter never drops a
  particle a brute-force scan says hits the wall; every kept lepton is above
  the analytic containment floor; the reduced parameter set runs.
- Units in `acc.dat`: `particles` in 1e10/bunch, `sigma_x/y` nm, `sigma_z` µm,
  `beta` mm, `emitt` 1e-6 m·rad. Internally positions are nm throughout.
