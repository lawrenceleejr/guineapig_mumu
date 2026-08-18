#!/usr/bin/env python3
"""Convert a GuineaPig bunch-crossing ("event") output into a HepMC2 ASCII
event file.

Inputs (all paths are optional except the main GuineaPig output file):
  --out       Main GuineaPig output file (e.g. muon_pairs_10tev_event1.out).
              Used to read the two beam energies (energy.1 / energy.2).
  --photons   photon.dat produced when store_photons/write_photons is set.
              Each line is "energy vx vy" with vx, vy the photon direction
              slopes (px/pz, py/pz) and a positive/negative energy sign used
              to indicate which beam the photon came from.
  --pairs     pairs.dat or pairs0.dat produced when track_pairs>0 or
              store_pairs>1. Each line is "energy vx vy vz [x y z]" with
              vx, vy, vz the normalised momentum direction (px/E, py/E, pz/E)
              of a pair-production lepton.
  --event-number   HepMC event number to record (default: 1).
  --output    Path to write the resulting HepMC2 ASCII file to.

The incoming muon beams are recorded as the two incoming particles of a
single vertex at the origin; every stored photon and pair-production lepton
is attached to that vertex as an outgoing (final-state) particle.  This is
not a full simulation of the underlying physics process, but it gives a
standard HepMC record of the particles GuineaPig produced for one bunch
crossing so it can be read by downstream HEP tools.
"""

import argparse
import math
import os
import re
import sys

MUON_MASS = 0.1056583715  # GeV
ELECTRON_MASS = 0.00051099895  # GeV
MUON_PDGID = 13
ELECTRON_PDGID = 11
PHOTON_PDGID = 22


def parse_beam_energies(out_path):
    """Read energy.1 / energy.2 (GeV) from a GuineaPig output file."""
    energy1 = energy2 = None
    if out_path and os.path.isfile(out_path):
        with open(out_path) as f:
            text = f.read()
        m1 = re.search(r"energy\.1\s*=\s*([-0-9.eE+]+)", text)
        m2 = re.search(r"energy\.2\s*=\s*([-0-9.eE+]+)", text)
        if m1:
            energy1 = float(m1.group(1))
        if m2:
            energy2 = float(m2.group(1))
    return energy1, energy2


def read_photons(path):
    """Yield (energy, vx, vy) tuples from a photon.dat file."""
    if not path or not os.path.isfile(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            yield float(parts[0]), float(parts[1]), float(parts[2])


def read_pairs(path):
    """Yield (energy, vx, vy, vz) tuples from a pairs.dat/pairs0.dat file."""
    if not path or not os.path.isfile(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            yield (float(parts[0]), float(parts[1]),
                   float(parts[2]), float(parts[3]))


def photon_four_vector(energy, vx, vy):
    e = abs(energy)
    pz = e / math.sqrt(1.0 + vx * vx + vy * vy)
    px = vx * pz
    py = vy * pz
    return px, py, pz, e


def pair_four_vector(energy, vx, vy, vz):
    e = abs(energy)
    px = vx * e
    py = vy * e
    pz = vz * e
    return px, py, pz, e


def build_particles(photon_path, pairs_path):
    """Return a list of (pdgid, px, py, pz, e, mass) outgoing particles."""
    particles = []
    for energy, vx, vy in read_photons(photon_path):
        px, py, pz, e = photon_four_vector(energy, vx, vy)
        particles.append((PHOTON_PDGID, px, py, pz, e, 0.0))

    for i, (energy, vx, vy, vz) in enumerate(read_pairs(pairs_path)):
        px, py, pz, e = pair_four_vector(energy, vx, vy, vz)
        # Pairs are produced as e+e-; alternate the charge assignment since
        # GuineaPig does not record it per-particle.
        pdgid = ELECTRON_PDGID if i % 2 == 0 else -ELECTRON_PDGID
        particles.append((pdgid, px, py, pz, e, ELECTRON_MASS))

    return particles


def write_hepmc(output_path, event_number, energy1, energy2, particles):
    energy1 = energy1 if energy1 is not None else 0.0
    energy2 = energy2 if energy2 is not None else 0.0

    beam1 = (MUON_PDGID, 0.0, 0.0, energy1, energy1, MUON_MASS)
    beam2 = (-MUON_PDGID, 0.0, 0.0, -energy2, energy2, MUON_MASS)

    n_out = len(particles)
    n_particles = 2 + n_out
    n_vertices = 1

    with open(output_path, "w") as f:
        f.write("HepMC::Version 2.06.09\n")
        f.write("HepMC::IO_GenEvent-START_EVENT_LISTING\n")
        f.write(
            "E %d -1 -1.0 -1.0 -1.0 0 0 %d 0 0 0 0\n"
            % (event_number, n_vertices)
        )
        f.write("U GEV MM\n")
        f.write("V -1 0 0.0 0.0 0.0 0.0 0 %d 0\n" % n_particles)

        barcode = 1
        for pdgid, px, py, pz, e, m in (beam1, beam2):
            f.write(
                "P %d %d %.9g %.9g %.9g %.9g %.9g 4 0.0 0.0 -1 0\n"
                % (barcode, pdgid, px, py, pz, e, m)
            )
            barcode += 1

        for pdgid, px, py, pz, e, m in particles:
            f.write(
                "P %d %d %.9g %.9g %.9g %.9g %.9g 1 0.0 0.0 0 0\n"
                % (barcode, pdgid, px, py, pz, e, m)
            )
            barcode += 1

        f.write("HepMC::IO_GenEvent-END_EVENT_LISTING\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", help="GuineaPig main output file")
    parser.add_argument("--photons", help="photon.dat file")
    parser.add_argument("--pairs", help="pairs.dat or pairs0.dat file")
    parser.add_argument("--event-number", type=int, default=1)
    parser.add_argument("--output", required=True, help="HepMC file to write")
    args = parser.parse_args()

    energy1, energy2 = parse_beam_energies(args.out)
    particles = build_particles(args.photons, args.pairs)
    write_hepmc(args.output, args.event_number, energy1, energy2, particles)

    print(
        "Wrote %d outgoing particle(s) to %s"
        % (len(particles), args.output)
    )


if __name__ == "__main__":
    sys.exit(main())
