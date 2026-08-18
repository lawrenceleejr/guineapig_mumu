# guineapig_mumu

This repository contains the **GuineaPig** beam-beam interaction simulation,
configured here for studies of a **10 TeV muon collider**.

GuineaPig was written by Daniel Schulte (DESY/CERN).  This README documents how
to build and run the code in this repository and summarises the input parameters
used for the reference 10 TeV muon pair-production test.

## Building

The default `makefile` points to a CERN AFS installation of FFTW that is not
available in most environments.  The CI and the recommended local build use the
`guinea_nofftw` target, which does not require FFTW:

```bash
make guinea_nofftw
```

Requirements: `gcc`, `make`, and the standard math library (`-lm`).

## Running

GuineaPig expects three command-line arguments:

```bash
./guinea_nofftw <accelerator> <parameter-set> <output-file>
```

The definitions for `<accelerator>` and `<parameter-set>` are read from
`acc.dat` (or the file `file_open` is hard-wired to read).  For example, the
10 TeV pair-production test is run as:

```bash
cat test_params.dat >> acc.dat
./guinea_nofftw muon muon_pairs_10tev muon_pairs_10tev.out
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
* `pairs.dat`   – stored pair-production particles (when `store_pairs=1` and
  `track_pairs=1`; not produced in the current no-tracking test).

## 10 TeV pair-production result

A reference run of `./guinea_nofftw muon muon_pairs_10tev muon_pairs_10tev.out`
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
and writes both the full run log (`guinea_pig.log`) and the GuineaPig output
file(s) to `./output` on the host, in addition to printing them to stdout.

To simulate more than one bunch crossing, pass the count as an argument or via
the `N_EVENTS` environment variable:

```bash
docker run --rm -v "$PWD/output":/output ghcr.io/lawrenceleejr/guineapig_mumu:latest 5
# or
docker run --rm -e N_EVENTS=5 -v "$PWD/output":/output ghcr.io/lawrenceleejr/guineapig_mumu:latest
```

The accelerator and parameter set can be overridden with the `ACCELERATOR` and
`PARAMS` environment variables (defaults: `muon` and `muon_pairs_10tev`).

To build the image locally:

```bash
docker build -t guineapig_mumu .
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
