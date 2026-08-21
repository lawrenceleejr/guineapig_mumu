#!/usr/bin/env python3
"""Convert a GuineaPig bunch-crossing ("event") output into a HepMC3 ASCII
(Asciiv3) event file.

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
  --output    Path to write the resulting HepMC3 ASCII file to.

The incoming muon beams are recorded as the two incoming particles of a
single vertex at the origin; every stored photon and pair-production lepton
is attached to that vertex as an outgoing (final-state) particle.  This is
not a full simulation of the underlying physics process, but it gives a
standard HepMC record of the particles GuineaPig produced for one bunch
crossing so it can be read by downstream HEP tools.

Reducing the record
-------------------
One bunch crossing is ~5.6e5 pair leptons, which is what makes the
downstream ddsim/MAIA job expensive.  The default filter here is geometric
and loss-free: a charged particle emitted at the IP into a solenoid field
follows a helix whose radial excursion from the axis is

    r(phi) = 2 * (pT / (0.3*B)) * sin(phi/2),    phi = 0.3*B*z/pz

(SI-ish units: metres, GeV, tesla).  If that never reaches the inner wall of
the beam-pipe/nozzle bore anywhere inside the detector, the particle stays in
vacuum from the IP to the end of the machine and can produce no hit and no
shower.  Such particles are dropped.

Because r(phi) saturates at 2*pT/(0.3*B) after a half turn, the crude version
of this test is the flat cut pT > 0.3*B*r_bore/2 (17.2 MeV for B=5 T and a
23 mm bore) -- that is what --pt-min does.  The flat cut is conservative but
wasteful: it keeps very forward, high-momentum pairs that only turn through a
few mrad inside the detector and so never get anywhere near the wall.  On a
measured bunch crossing, pT > 15 MeV keeps 238038 leptons; the helix test
against a 23 mm + 10 mm/m bore keeps 202537 as implemented here, or 196048 if
the off-axis start is treated exactly instead of conservatively (see
reaches_bore).  The exact set is a strict subset of what the flat cut keeps --
it is the same set minus 41990 particles that provably never leave the bore --
while the conservative one additionally keeps 582 leptons that start far
enough off-axis to need less than 15 MeV to reach the wall.

--bore-r0 / --bore-slope / --bore-zmax describe that bore and MUST be set from
the geometry actually being simulated (the beam-pipe and nozzle solids in the
MAIA compact XML).  The defaults describe a 23 mm bore opening at 10 mm/m out
to |z| = 6 m, which reproduces the 17 MeV figure quoted for MAIA but has not
been checked against the XML.  A bore narrower than the value used here would
make the filter lossy.
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

# Which sign of the energy field in pairs0.dat/pairs.dat corresponds to the
# electron.  GuineaPig guarantees the two members of a pair have opposite
# signs but the absolute species assignment is a convention; flip this single
# flag if the downstream convention differs.
POSITIVE_ENERGY_IS_ELECTRON = True

# Default transverse-momentum cut applied to the charged pair-production
# leptons, in GeV.  GuineaPig itself applies no pt cut: store_pair() cuts on
# ENERGY (pair_ecut, default 5 MeV) and the pt>20 MeV / theta>0.15 test in
# background.c only feeds the printed anzahl_* counters, never the stored
# output.  So without this the HepMC keeps pairs down to ~1 keV of pt.
#
# For reference, the MAIA inside-beam-pipe value is 17 MeV: below that a
# charged particle's helix in the solenoid field keeps it inside the beam
# pipe, so it never reaches the detector.  The default here is 15 MeV,
# slightly looser so that the beam-pipe threshold can still be applied
# downstream without having to regenerate the sample.
#
# How much this cut costs depends entirely on WHICH pair file it is applied to:
#   pairs0.dat (production time) : pT>15 MeV keeps ~1% of pairs
#   pairs.dat  (tracked)         : pT>15 MeV keeps ~41% of pairs
# On the tracked sample this is therefore a consequential physics choice, not a
# cheap cleanup -- it removes most of the pairs that actually reach the detector.
# It is also NOT a useful way to shrink the downstream ddsim output: that file is
# dominated by the MCParticle collection (>99% of bytes for a cut sample), not by
# hits, so SIM.part.minimalKineticEnergy is the lever for file size.  Pass 0 to
# disable this cut.
DEFAULT_PT_MIN = 0.0

# Geometry of the vacuum channel the pairs are born into, used by the
# loss-free helix filter (see the module docstring).  These describe the inner
# wall of the beam-pipe/nozzle bore as r_bore(z) = BORE_R0 + BORE_SLOPE*|z|,
# valid out to |z| = BORE_ZMAX.  They are NOT taken from the MAIA XML -- set
# them from the geometry you actually simulate.
DEFAULT_B_FIELD = 5.0        # tesla, MAIA solenoid
DEFAULT_BORE_R0 = 23.0       # mm, bore radius at the IP
DEFAULT_BORE_SLOPE = 10.0    # mm per metre of |z|
DEFAULT_BORE_ZMAX = 6.0      # m, end of the nozzle / simulated region


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
    """Yield (energy, vx, vy, vz, x, y) tuples from a pairs.dat/pairs0.dat file.

    Columns 5-7 are the particle position; for pairs.dat that is the position
    at the end of the beam-beam tracking, in GuineaPig's internal transverse
    unit of nm (cut_z is printed in um but held internally in nm too, see the
    grid->cut_z*1e-3 in guinea_pig.c).  x and y are returned in metres; z is
    ignored because it never exceeds ~6 mm, which is negligible against the
    metre-scale bore length.  Files without the position columns yield 0.0.
    """
    if not path or not os.path.isfile(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            x = float(parts[4]) * 1e-9 if len(parts) >= 6 else 0.0
            y = float(parts[5]) * 1e-9 if len(parts) >= 6 else 0.0
            yield (float(parts[0]), float(parts[1]),
                   float(parts[2]), float(parts[3]), x, y)


def reaches_bore(pt, pz, r_prod, b_field, bore_r0, bore_slope, bore_zmax):
    """True if a charged particle's helix reaches the bore wall in the detector.

    Solves, in closed form, whether there is a turn angle phi in (0, phi_max]
    with

        2*R*sin(phi/2) >= r_bore(z(phi)) - r_prod,
        R = pt/(0.3*B),  z(phi) = |pz|*phi/(0.3*B),
        r_bore(z) = bore_r0 + bore_slope*z

    i.e. whether the chord from the production point ever spans the gap to the
    wall.  Subtracting the production radius r_prod rather than doing the
    trigonometry with a real off-axis start is deliberately conservative: it
    assumes the offset points straight at the wall, so the filter can only
    ever keep a particle it should have kept, never drop one.  Measured on one
    bunch crossing against an exact off-axis trajectory scan: the exact test
    keeps 196048 leptons, this one keeps 202537, so the safety costs 3.3%
    extra particles.  Ignoring r_prod altogether would keep 195552, i.e. would
    silently lose 496 particles that do reach the wall -- which is the thing
    this filter exists not to do.  The exact test is not used because r(phi)
    with an off-axis start is not concave, so it needs a per-particle scan,
    and this script has to run in the runtime image without numpy.

    Units: metres, GeV, tesla, radians.

    f(phi) = 2*R*sin(phi/2) - (bore_r0 - r_prod) - bore_slope*|pz|/(0.3*B)*phi
    is concave on [0, pi] and strictly decreasing beyond it (the sine turns
    over while the wall keeps receding), so its maximum over all phi > 0 lies
    at the stationary point phi* = 2*acos(min(1, slope*|pz|/pt)), clamped to
    the range of phi actually spent inside the simulated region.  Evaluating f
    once at that point therefore decides the question exactly -- no scan.
    """
    gap = bore_r0 - r_prod                      # wall distance to span, m
    if gap <= 0.0:
        return True                             # born on/outside the wall
    if pt <= 0.0:
        return False                            # never leaves the axis
    omega = 0.3 * b_field                       # rad per metre, per GeV of pz
    chord_max = 2.0 * pt / omega                # max radial excursion, m
    apz = abs(pz)
    # phi available before the particle leaves the simulated region
    phi_lim = math.pi if apz <= 0.0 else min(math.pi, omega * bore_zmax / apz)
    # rate at which the wall recedes, per radian of turn
    recede = bore_slope * apz / omega
    if recede <= 0.0:
        phi_star = math.pi
    else:
        ratio = 2.0 * recede / chord_max        # = bore_slope*|pz|/pt
        if ratio >= 1.0:
            return False                        # wall outruns the helix at once
        phi_star = 2.0 * math.acos(ratio)
    phi = min(phi_star, phi_lim)
    return chord_max * math.sin(0.5 * phi) >= gap + recede * phi


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


def build_particles(photon_path, pairs_path, pt_min=DEFAULT_PT_MIN,
                    b_field=DEFAULT_B_FIELD, bore_r0=DEFAULT_BORE_R0,
                    bore_slope=DEFAULT_BORE_SLOPE,
                    bore_zmax=DEFAULT_BORE_ZMAX, geometric=True,
                    keep_photons=True):
    """Return (particles, counts).

    particles is a list of (pdgid, px, py, pz, e, mass) outgoing particles.
    counts is a dict of bookkeeping numbers for the log line.

    Both cuts apply only to the charged pair leptons.  The solenoid does not
    bend photons, so a photon that starts inside the bore stays inside it and
    the geometric argument would drop every one of them; they are governed by
    --keep-photons instead.
    """
    particles = []
    counts = {"photons": 0, "pairs_in": 0, "cut_pt": 0, "cut_geom": 0}

    if keep_photons:
        for energy, vx, vy in read_photons(photon_path):
            px, py, pz, e = photon_four_vector(energy, vx, vy)
            particles.append((PHOTON_PDGID, px, py, pz, e, 0.0))
            counts["photons"] += 1

    # mm -> m for the bore description
    r0_m = bore_r0 * 1e-3
    slope_m = bore_slope * 1e-3

    for energy, vx, vy, vz, x, y in read_pairs(pairs_path):
        px, py, pz, e = pair_four_vector(energy, vx, vy, vz)
        counts["pairs_in"] += 1
        pt = math.hypot(px, py)
        if pt < pt_min:
            counts["cut_pt"] += 1
            continue
        if geometric and not reaches_bore(pt, pz, math.hypot(x, y), b_field,
                                          r0_m, slope_m, bore_zmax):
            counts["cut_geom"] += 1
            continue
        # GuineaPig records the charge in the SIGN of the energy: store_full_pair()
        # in background.c negates the second member of every pair ("e2 = -e2"), so
        # the two members always carry opposite signs.  Index parity is NOT usable:
        # pair_ecut and pair_ratio drop members individually inside store_pair(),
        # which desynchronises the alternation (measured agreement with the true
        # sign on a real pairs0.dat: 49.5%, i.e. a coin flip).
        pdgid = ELECTRON_PDGID if (energy > 0.0) == POSITIVE_ENERGY_IS_ELECTRON \
            else -ELECTRON_PDGID
        particles.append((pdgid, px, py, pz, e, ELECTRON_MASS))

    return particles, counts


def write_hepmc(output_path, event_number, energy1, energy2, particles,
                weight=1.0):
    """Write the HepMC3 ASCII record.

    weight is the HepMC event weight and is how a GuineaPig pair_ratio < 1 run
    stays interpretable: with pair_ratio=r each stored lepton stands for 1/r
    real ones, so downstream hit densities must be multiplied by 1/r.  ddsim
    does NOT act on this -- it is bookkeeping for the analysis step, carried
    here so the thinning factor cannot get lost between the two.
    """
    energy1 = energy1 if energy1 is not None else 0.0
    energy2 = energy2 if energy2 is not None else 0.0

    beam1 = (MUON_PDGID, 0.0, 0.0, energy1, energy1, MUON_MASS)
    beam2 = (-MUON_PDGID, 0.0, 0.0, -energy2, energy2, MUON_MASS)

    n_out = len(particles)
    n_particles = 2 + n_out
    n_vertices = 1
    vertex_id = -1

    with open(output_path, "w") as f:
        f.write("HepMC::Version 3.02.05\n")
        f.write("HepMC::Asciiv3-START_EVENT_LISTING\n")
        f.write("E %d %d %d\n" % (event_number, n_vertices, n_particles))
        f.write("U GEV MM\n")
        f.write("W %.9g\n" % weight)

        barcode = 1
        beam_ids = []
        for pdgid, px, py, pz, e, m in (beam1, beam2):
            f.write(
                "P %d 0 %d %.9g %.9g %.9g %.9g %.9g 4\n"
                % (barcode, pdgid, px, py, pz, e, m)
            )
            beam_ids.append(barcode)
            barcode += 1

        f.write(
            "V %d 0 [%s]\n" % (vertex_id, ",".join(str(i) for i in beam_ids))
        )

        for pdgid, px, py, pz, e, m in particles:
            f.write(
                "P %d %d %d %.9g %.9g %.9g %.9g %.9g 1\n"
                % (barcode, vertex_id, pdgid, px, py, pz, e, m)
            )
            barcode += 1

        f.write("HepMC::Asciiv3-END_EVENT_LISTING\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", help="GuineaPig main output file")
    parser.add_argument("--photons", help="photon.dat file")
    parser.add_argument("--pairs", help="pairs.dat or pairs0.dat file")
    parser.add_argument("--event-number", type=int, default=1)
    parser.add_argument(
        "--pt-min", type=float, default=DEFAULT_PT_MIN, metavar="GEV",
        help="flat transverse-momentum cut in GeV on the charged pair leptons, "
             "applied in addition to the geometric filter (default: "
             "%(default)s, i.e. off). This is the crude form of the same "
             "containment argument -- 0.3*B*r_bore/2, or 0.017 for a 23 mm "
             "bore at 5 T -- and is kept only for reproducing older samples.")
    parser.add_argument(
        "--no-geometric-filter", dest="geometric", action="store_false",
        help="keep every pair lepton regardless of whether its helix can ever "
             "reach the bore wall. Makes the file ~2.9x bigger.")
    parser.add_argument(
        "--b-field", type=float, default=DEFAULT_B_FIELD, metavar="TESLA",
        help="solenoid field used by the geometric filter (default: "
             "%(default)s, the MAIA value).")
    parser.add_argument(
        "--bore-r0", type=float, default=DEFAULT_BORE_R0, metavar="MM",
        help="inner radius of the beam-pipe/nozzle bore at the IP (default: "
             "%(default)s). SET THIS FROM THE GEOMETRY YOU SIMULATE -- a real "
             "bore narrower than this makes the filter lossy.")
    parser.add_argument(
        "--bore-slope", type=float, default=DEFAULT_BORE_SLOPE,
        metavar="MM_PER_M",
        help="rate at which the bore opens with |z| (default: %(default)s). "
             "0 treats it as a cylinder, which is the conservative choice.")
    parser.add_argument(
        "--bore-zmax", type=float, default=DEFAULT_BORE_ZMAX, metavar="M",
        help="|z| at which the particle has left the simulated region "
             "(default: %(default)s).")
    parser.add_argument(
        "--drop-photons", dest="keep_photons", action="store_false",
        help="omit the beamstrahlung photons. They are 0.3%% of the record and "
             "every one measured lies within 3 mrad of the axis, so none can "
             "reach a detector element; they are also not weighted up to "
             "physical rates the way the pairs are (store_photon() carries no "
             "macroparticle weight, store_full_pair() does), so they should "
             "not be used for a photon-background estimate.")
    parser.add_argument(
        "--weight", type=float, default=1.0,
        help="HepMC event weight. Set to 1/pair_ratio when GuineaPig was run "
             "with pair_ratio < 1 (default: %(default)s).")
    parser.add_argument("--output", required=True, help="HepMC file to write")
    args = parser.parse_args()

    energy1, energy2 = parse_beam_energies(args.out)
    particles, counts = build_particles(
        args.photons, args.pairs, args.pt_min, args.b_field, args.bore_r0,
        args.bore_slope, args.bore_zmax, args.geometric, args.keep_photons)
    write_hepmc(args.output, args.event_number, energy1, energy2, particles,
                args.weight)

    kept_pairs = len(particles) - counts["photons"]
    print(
        "Wrote %d outgoing particle(s) to %s: %d photon(s) + %d/%d pair "
        "lepton(s)" % (len(particles), args.output, counts["photons"],
                       kept_pairs, counts["pairs_in"])
    )
    if counts["cut_pt"]:
        print("  flat pT > %g GeV dropped %d" % (args.pt_min, counts["cut_pt"]))
    if args.geometric:
        print("  geometric filter dropped %d (B=%g T, bore = %g mm + %g mm/m "
              "* |z| out to |z| < %g m; these never leave the bore)"
              % (counts["cut_geom"], args.b_field, args.bore_r0,
                 args.bore_slope, args.bore_zmax))
    if args.weight != 1.0:
        print("  event weight %g -- multiply downstream hit densities by this"
              % args.weight)


if __name__ == "__main__":
    sys.exit(main())
