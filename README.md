# guineapig_mumu

This repository contains the **GuineaPig** beam-beam interaction simulation,
configured here for studies of a **10 TeV muon collider**.

GuineaPig was written by Daniel Schulte (DESY/CERN).  This README documents how
to build and run the code in this repository and summarises the input parameters
used for the reference 10 TeV muon pair-production test.

## Building

Build the `guinea` target, which uses the FFT Poisson solver (`fourtrans3.c`):

```bash
sudo apt-get install -y gcc make libfftw3-dev   # Debian/Ubuntu
make guinea FFTW_HOME=/usr
```

`FFTW_HOME` defaults to `/usr` and can be pointed anywhere (it used to be
hard-wired to a CERN AFS path that no longer resolves).  The link line sets
`-rpath`, so the binary finds `libfftw3` without `LD_LIBRARY_PATH`.

A `guinea_nofftw` target still exists and needs no FFTW, but it omits
`fourtrans3.c`, so the two builds compute the beam-beam field differently.
**Do not mix results from the two builds in one comparison.**

Requirements: `gcc`, `make`, `libfftw3-dev`, and the standard math library.

## Running

GuineaPig expects three command-line arguments:

```bash
./guinea <accelerator> <parameter-set> <output-file>
```

The definitions for `<accelerator>` and `<parameter-set>` are read from
`acc.dat` (or the file `file_open` is hard-wired to read).  For example, the
10 TeV pair-production test is run as:

```bash
cat test_params.dat >> acc.dat
./guinea muon muon_pairs_10tev muon_pairs_10tev.out
```

The `cat test_params.dat >> acc.dat` step appends the CI test parameter blocks
to the main input file.

## Units

All accelerator variables in `acc.dat` use the following units (see the
GuineaPig manual):

| Variable    | Unit                                   |
|-------------|----------------------------------------|
| `particles` | number of particles per bunch in **1e10** |
| `energy`    | GeV                                    |
| `emitt_x/y` | normalised emittance in 1e-6 m rad     |
| `beta_x/y`  | mm                                     |
| `sigma_x/y` | nm                                     |
| `sigma_z`   | um                                     |
| `offset_x/y`| nm                                     |
| `waist_x/y` | um                                     |
| `angle_x/y` | rad                                    |
| `charge_sign`| -1 for opposite-sign beams, +1 for same-sign |

**Important:** `particles` is in units of 1e10 particles per bunch.  The default
`muon` accelerator has `particles=180.0`, i.e. `1.8e12` particles per bunch.

## Accelerator definitions (`acc.dat`)

Three reference muon-collider accelerators are provided:

* `muon`  – baseline 10 TeV configuration: 5000 GeV/beam, 1.8e12 particles/bunch,
  beta*=1.5 mm, normalized emittance 25e-6 m rad.
* `muon2` – 5000 GeV/beam with 0.707 of the baseline bunch population.
* `muon4` – 5000 GeV/beam with 0.5 of the baseline bunch population.

## Parameter sets

The repository ships two small parameter sets in `test_params.dat` that are used
by the CI workflow:

* `ci_test` – minimal fast luminosity-only smoke test.
* `muon_pairs_10tev` – 10 TeV muon-collider pair-production test.

### `muon_pairs_10tev` settings

| Parameter       | Value | Meaning |
|-----------------|-------|---------|
| `n_x`, `n_y`    | 32    | transverse grid cells |
| `n_z`           | 20    | longitudinal slices |
| `n_t`           | 3     | timesteps per slice crossing |
| `n_m`           | 5000  | macroparticles per beam |
| `do_eloss`      | 1     | include beamstrahlung energy loss |
| `do_photons`    | 1     | generate/store beamstrahlung photons |
| `do_pairs`      | 1     | enable pair production |
| `store_pairs`   | 1     | store pair production information |
| `track_pairs`   | 0     | do not track pairs through fields |
| `grids`         | 1     | number of grids for tracking (1 because tracking is off) |
| `ecm_min_gg`    | 9990.0 | photon-photon threshold near 10 TeV |
| `beam_size`     | 0     | do not include beam-size effect on luminosity |
| `silent`        | 1     | reduce printed output |

`track_pairs=0` is used because the `guinea_nofftw` build crashes when tracking
pairs through multiple grids.  Pair production is still evaluated and counted.

## Outputs

The main output file (`<output-file>`) contains the beam, grid, and result
parameters.  Additional files may be produced depending on the switches:

* `lumi.ee.out` – energies for e+e- scatters (when `do_lumi=1`).
* `photon.dat`  – stored beamstrahlung photons (when `do_photons=1`).
* `pairs.dat`/`pairs0.dat` – stored pair-production particles (when
  `store_pairs>=1`; not produced in the current no-tracking test since
  `track_pairs=0`).
* `<params>_event<n>.hepmc` – a HepMC3 ASCII event file for each bunch
  crossing, produced automatically as a post-processing step by
  `docker/make_hepmc.py`. It records the two incoming beams (read from
  `energy.1`/`energy.2` in the main output file) plus every photon in
  `photon.dat` and, when available, every pair-production lepton in
  `pairs.dat`/`pairs0.dat` as final-state particles of a single vertex.
* `<params>.hepmc` – the per-event HepMC3 files for the whole run merged into
  a single multi-event HepMC3 file, produced automatically by
  `docker/merge_hepmc.py` (see [Docker](#docker) below).

## 10 TeV pair-production result

A reference run of `./guinea muon muon_pairs_10tev muon_pairs_10tev.out`
is committed as `muon_pairs_10tev.out`.  For this single bunch crossing (one
event) with 1.8e12 muons per bunch, the output reports:

```
pairs 943808 692645
< E_cm > 10000 0.668872
```

The two numbers after `pairs` are the numbers of produced pairs recorded by the
two internal bookkeeping counters; the `< E_cm >` line confirms the
centre-of-mass energy is 10 TeV.

## Docker

A Docker image with `guinea_nofftw` and its input files pre-built is published
automatically by CI to the GitHub Container Registry at
`ghcr.io/lawrenceleejr/guineapig_mumu`. No setup is required beyond `docker run`:

```bash
docker run --rm -v "$PWD/output":/output ghcr.io/lawrenceleejr/guineapig_mumu:latest
```

This simulates one bunch crossing ("event") of the `muon_pairs_10tev` process
and writes the full run log (`guinea_pig.log`), the GuineaPig output
file(s), and a HepMC3 event file (`<params>_event<n>.hepmc`) to `./output` on
the host, in addition to printing them to stdout. No extra flags are needed
to get the HepMC file; it is produced automatically for every run.

To simulate more than one bunch crossing, pass the count as an argument or via
the `N_EVENTS` environment variable:

```bash
docker run --rm -v "$PWD/output":/output ghcr.io/lawrenceleejr/guineapig_mumu:latest 5
# or
docker run --rm -e N_EVENTS=5 -v "$PWD/output":/output ghcr.io/lawrenceleejr/guineapig_mumu:latest
```

By default, the per-event HepMC files for the run are also merged into a
single `<params>.hepmc` file (via `docker/merge_hepmc.py`), so a `N_EVENTS=5`
run produces `muon_pairs_10tev_event1.hepmc` .. `muon_pairs_10tev_event5.hepmc`
as well as a combined `muon_pairs_10tev.hepmc` containing all 5 events. Set
`MERGE_HEPMC=0` to skip this and keep only the per-event files.

The accelerator and parameter set can be overridden with the `ACCELERATOR` and
`PARAMS` environment variables (defaults: `muon` and `muon_pairs_10tev`).

To build the image locally:

```bash
docker build -t guineapig_mumu .
```

### Singularity / Apptainer

The same image can be run with [Singularity](https://sylabs.io/singularity/)
or its drop-in replacement [Apptainer](https://apptainer.org/) (the
`singularity` command below works identically with `apptainer`), pulling
directly from GHCR:

```bash
mkdir -p output
singularity run --bind "$PWD/output":/output docker://ghcr.io/lawrenceleejr/guineapig_mumu:latest
```

Environment variables are passed through with `SINGULARITYENV_` (or
`APPTAINERENV_`) prefixes, and positional arguments are given after `--`:

```bash
singularity run --bind "$PWD/output":/output \
    --env ACCELERATOR=muon,PARAMS=muon_pairs_10tev,N_EVENTS=5 \
    docker://ghcr.io/lawrenceleejr/guineapig_mumu:latest
```

To build a local `.sif` image file instead of pulling on every run:

```bash
singularity build guineapig_mumu.sif docker://ghcr.io/lawrenceleejr/guineapig_mumu:latest
singularity run --bind "$PWD/output":/output guineapig_mumu.sif
```

## CI

`.github/workflows/ci.yml` builds `guinea_nofftw` on every push and pull request
and runs the `ci_test` and `muon_pairs_10tev` parameter sets. On pushes to
`main` and version tags (`v*`), it also builds the Docker image and publishes
it to GHCR.

## Reference

GuineaPig manual (SLAC NLC webpage version), K. Thompson revision of the
original manual by D. Schulte:
https://indico.ihep.ac.cn/event/26852/contributions/196716/attachments/93232/122059/GuineaPigManual.pdf


## Pair output: use `pairs.dat`, not `pairs0.dat`

GuineaPig writes two different pair files and they are **not** interchangeable.

| file | written by | contents |
|---|---|---|
| `pairs0.dat` | `store_pair()`, gated on `store_pairs>1` | pair kinematics **at production**, before the opposing beam's field deflects them |
| `pairs.dat` | `print_pairs()`, requires `track_pairs>0` | pair kinematics **after** tracking through the beam-beam field |

The beam-field deflection is what gives incoherent pairs their transverse
momentum, so `pairs0.dat` is not usable for detector background studies.
Measured on one bunch crossing (560,566 pair leptons, identical particles in
both files):

| | production (`pairs0.dat`) | tracked (`pairs.dat`) | ratio |
|---|---|---|---|
| mean pT | 1.730 MeV | 39.063 MeV | x22.6 |
| median pT | 0.665 MeV | 9.833 MeV | x14.8 |
| N(pT > 18 MeV) | 5,811 (1.04%) | 219,510 (39.2%) | **x37.8** |
| mean energy | 1.305 GeV | 1.141 GeV | x0.874 (radiated) |

18 MeV is where an on-axis particle's radial excursion `2*pT/(0.3*B)` reaches the
24 mm vertex layer in a 5 T solenoid.  Untracked pairs stay inside the beam pipe
(median excursion 0.89 mm) and produce essentially no detector hits; tracked ones
reach 13 mm median.

`muon_pairs_10tev` therefore sets `track_pairs=1`.  The previous configuration is
kept as `muon_pairs_10tev_notrack` **only** for reproducing old samples.

### Two fixes were needed to make tracking work

1. **`step_pair_1()` buffer overflow.** `synrad()` emits up to 10000 photons per
   call (its own guard is `j>=10000`) and the two beam-tracking callers declare
   `photon[10000]`, but the pair-tracking caller declared `ph[1000]`.  A pair
   radiating >1000 photons in one step wrote past the end of that stack array,
   aborting with `*** stack smashing detected ***`.  The adjacent `vx0/vy0/vz0`
   were corrupted first, which surfaced as `|v|>1` in the energy-consistency
   check.  Fixed to `ph[10000]`.
2. **`pair_step=5.0`** (not the `0.2` default).  `d_eps` is proportional to
   `pair_step` and the sub-step is `step0/(d+1)`, so a *larger* `pair_step`
   subdivides more finely.  At `0.2` a pair exceeds synrad's 10000-photon limit
   in one step and the run aborts with "too many photons produced by one
   particle".  **This value has not been checked for convergence** -- vary it
   before trusting absolute numbers.

The CI asserts that `pairs.dat` exists and that its mean pT is well above
`pairs0.dat`'s, so this cannot silently regress.

## Docker

The image builds the FFTW `guinea` target and runs the tracked configuration:

```bash
docker run --rm -v "$PWD/out:/output" ghcr.io/lawrenceleejr/guineapig_mumu 1
```

Environment variables: `ACCELERATOR` (default `muon`), `PARAMS`
(`muon_pairs_10tev`), `N_EVENTS` (1), `PT_MIN` (0.015), `OUTPUT_DIR`
(`/output`).

If `pairs.dat` is missing the entrypoint now prints a loud warning before
falling back to `pairs0.dat`; the silent fallback is how the untracked samples
went unnoticed.

## Making the sample cheap enough to simulate in MAIA

One bunch crossing is **560,566 pair leptons and 1,810 photons**.  The leptons
are 99.7% of the record and effectively all of the downstream `ddsim` cost.
They are soft and very forward: median energy 68 MeV, median pT 9.8 MeV,
median polar angle 6.6 deg.

Three independent levers, in decreasing order of how safe they are.  All
numbers below are measured on one tracked bunch crossing of
`muon_pairs_10tev`.

### 1. Geometric filter -- loss-free, on by default (x2.77)

A charged particle leaving the IP into the solenoid follows a helix whose
distance from the beam axis is

```
r(phi) = 2 * (pT / (0.3*B)) * sin(phi/2),     phi = 0.3*B*z/pz
```

(metres, GeV, tesla).  If that never reaches the inner wall of the
beam-pipe/nozzle bore anywhere inside the simulated region, the particle flies
from the IP to the end of the machine in vacuum: no hit, and no shower in the
tungsten either.  Dropping it changes nothing downstream.

`docker/make_hepmc.py` applies this test in closed form.  Because `r(phi)` is
concave in `phi` up to a half turn and strictly falling after it, the maximum
excursion is at a single stationary point, so one evaluation decides each
particle -- 560k leptons in 2 s with no numpy.

| filter | leptons kept | HepMC size |
|---|---|---|
| none | 560,566 | 44.0 MB |
| flat `pT > 15 MeV` (the old default) | 238,038 | 18.1 MB |
| geometric, 23 mm + 10 mm/m bore | 202,537 | 15.5 MB |

The flat pT cut is the crude form of the same argument: `r(phi)` saturates at
`2*pT/(0.3*B)`, so requiring the particle to be *able* to reach `r_bore` gives
`pT > 0.3*B*r_bore/2` = 17.2 MeV for a 23 mm bore at 5 T.  That is where the
old `PT_MIN` came from.  It is conservative but wasteful, because it ignores
*where* the particle gets there: 41,990 of the leptons it keeps are so forward
that they turn through only a few mrad before leaving the detector and never
come near the wall.  Those carry most of the pair energy -- the filter keeps
38% of the energy where the flat cut keeps 86% -- because the high-energy pairs
are precisely the collinear ones that go straight down the bore.  They are the
machine's problem in the MDI, not a source of detector hits.

**`BORE_R0` / `BORE_SLOPE` / `BORE_ZMAX` must be set from the geometry you
actually simulate** (the beam-pipe and nozzle solids in the MAIA compact XML).
The defaults -- 23 mm at the IP, opening 10 mm/m, out to |z| = 6 m -- reproduce
the 17 MeV figure quoted for MAIA and are consistent with the vertex barrel
starting at r = 3.0 cm and the nozzle shadowing |eta| > 2.44, but they have not
been checked against the XML.  A real bore *narrower* than `BORE_R0` makes the
filter lossy, so err small.  Sensitivity, at `BORE_SLOPE=10`:

| `BORE_R0` | leptons kept |
|---|---|
| 10 mm | 246,276 |
| 15 mm | 223,760 |
| 23 mm | 196,048 |
| 30 mm | 175,998 |

(Those are with an exact off-axis trajectory; the shipped filter is
deliberately a little looser, see below.)

The filter treats the production offset conservatively: it requires the helix
chord to span `r_bore - r_prod` rather than doing the trigonometry with a real
off-axis start, so it can only ever keep a particle it should have kept.  On
this crossing the exact test keeps 196,048 and the shipped one 202,537, i.e.
the safety costs 3.3% extra particles.  Ignoring the offset entirely would keep
195,552 -- it would silently lose 496 particles that do reach the wall.

### 2. Cut at the source with `pair_ecut` (x1.35, ~0.7% loss)

`pair_ecut` is GuineaPig's own energy cut, applied in `store_pair()` before the
particle joins the tracking list.  A lepton needs at least 17.2 MeV of *energy*
before it can have 17.2 MeV of pT, so raising `pair_ecut` from the 5 MeV
default to 0.017 is nearly free: the stored sample drops from 560,566 to
415,515 while the number surviving the geometric filter changes from 196,048 to
194,743, i.e. 0.7%.  It is not exactly free because the cut is applied at
production, before the opposing beam's field has done work on the particle.

### 3. Thin with `pair_ratio` (arbitrary factor, unbiased, costs statistics)

`pair_ratio=0.1` stores a random 10% of the pairs.  It is applied in
`store_pair()` and thins uniformly -- scaled by 10, the sample reproduces the
full one to 0.1% in count and 1.4% in energy:

| | leptons stored | after geometric filter, x1/ratio |
|---|---|---|
| `pair_ratio=1` | 560,566 | 196,048 |
| `pair_ratio=0.1` | 56,483 | 196,220 |

Every downstream hit count must then be multiplied by `1/pair_ratio`.
`make_hepmc.py` writes that as the HepMC event weight (`W` record) and the
Docker entrypoint derives it from the `pair_ratio` GuineaPig reports, so it
cannot get lost between the two -- but nothing applies it for you, and `ddsim`
ignores it.

**This is the one lever that is not safe by default.**  Scaling up 10% of the
pairs gives the right *mean* occupancy with the wrong fluctuations, so it is
fine for hit-density maps and wrong for anything that depends on the true local
density in a single crossing: cluster merging, tracking confusion, occupancy
tails.  For those, use the geometric filter alone and push the remaining cost
into `ddsim` settings.

### Stacking them

`muon_pairs_10tev_reduced` is levers 2 and 3 together:

| configuration | leptons in the HepMC | size |
|---|---|---|
| `muon_pairs_10tev`, no filter | 560,566 | 44.0 MB |
| `muon_pairs_10tev`, old `PT_MIN=0.015` | 238,038 | 18.1 MB |
| `muon_pairs_10tev`, geometric filter | 202,537 | 15.5 MB |
| `muon_pairs_10tev_reduced`, geometric filter, `W=10` | 19,549 | 1.5 MB |

```bash
docker run --rm -v "$PWD/out:/output" \
    -e PARAMS=muon_pairs_10tev_reduced \
    -e BORE_R0=23.0 -e BORE_SLOPE=10.0 \
    ghcr.io/lawrenceleejr/guineapig_mumu 1
```

### What this does not fix

* **`ddsim` CPU is spent on showers, not on primaries.**  The leptons the filter
  keeps are exactly the ones that hit tungsten, so the remaining cost is the
  nozzle shower.  `SIM.part.minimalKineticEnergy` and region-based range cuts in
  the nozzle are the levers there, and they are independent of anything in this
  repository.  If the same background is being re-simulated for every physics
  event, simulating it once into a reusable overlay library will beat any cut
  here.
* **The photons are not usable as a photon background.**  All 1,810 lie within
  3 mrad of the axis, so none can reach a detector element, and they are not
  normalised the way the pairs are: `store_full_pair()` samples against the
  macroparticle weight `wgt` so the pairs come out at physical rates, while
  `store_photon()` carries no weight at all.  `DROP_PHOTONS=1` removes them;
  they are only 0.3% of the record either way.
* **Production vertices are discarded.**  Every particle is written at the
  origin, though `pairs.dat` records where the beam-beam tracking left it
  (columns 5-7, in nm: median 0.18 mm from the axis, max 5.6 mm).  That is a
  0.25% effect on the filter, but it does smear the background's true d0/z0
  distribution in the vertex detector.
