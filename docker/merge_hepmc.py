#!/usr/bin/env python3
"""Merge per-event HepMC2 ASCII files into a single multi-event file.

make_hepmc.py writes one file per bunch crossing.  A HepMC2 IO_GenEvent
stream can hold any number of events between a single START/END listing
pair, so merging is a matter of keeping one header, concatenating the event
bodies, and emitting one END marker:

    HepMC::Version 2.06.09
    HepMC::IO_GenEvent-START_EVENT_LISTING
    E 1 ...  <- event 1 block (E, U, V, P... lines)
    E 2 ...  <- event 2 block
    ...
    HepMC::IO_GenEvent-END_EVENT_LISTING

Inputs are sorted by the event number embedded in the filename, so
event10 sorts after event9 rather than after event1.

Usage:
    merge_hepmc.py --output all.hepmc <files...>
    merge_hepmc.py --output all.hepmc --renumber <files...>
"""

import argparse
import os
import re
import sys

VERSION_PREFIX = "HepMC::Version"
START_MARKER = "HepMC::IO_GenEvent-START_EVENT_LISTING"
END_MARKER = "HepMC::IO_GenEvent-END_EVENT_LISTING"


def event_sort_key(path):
    """Sort by the integer in 'event<N>' if present, else by name."""
    m = re.search(r"event(\d+)", os.path.basename(path))
    return (0, int(m.group(1))) if m else (1, os.path.basename(path))


def merge(paths, output_path, renumber=False):
    n_events = 0
    n_particles = 0
    version_line = "%s 2.06.09\n" % VERSION_PREFIX

    with open(output_path, "w") as out:
        out.write(version_line)
        out.write(START_MARKER + "\n")

        for path in paths:
            with open(path) as f:
                for line in f:
                    if line.startswith(VERSION_PREFIX) or \
                            line.startswith(START_MARKER) or \
                            line.startswith(END_MARKER):
                        continue
                    if line.startswith("E "):
                        n_events += 1
                        if renumber:
                            parts = line.split(" ", 2)
                            line = "E %d %s" % (n_events, parts[2])
                    elif line.startswith("P "):
                        n_particles += 1
                    out.write(line)

        out.write(END_MARKER + "\n")

    return n_events, n_particles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="per-event HepMC files")
    parser.add_argument("--output", required=True, help="merged file to write")
    parser.add_argument(
        "--renumber", action="store_true",
        help="renumber events sequentially from 1 (default: keep the event "
             "numbers already recorded in each file)")
    args = parser.parse_args()

    paths = sorted(args.inputs, key=event_sort_key)
    n_events, n_particles = merge(paths, args.output, args.renumber)
    print("Merged %d file(s) -> %s: %d event(s), %d particle(s)"
          % (len(paths), args.output, n_events, n_particles))


if __name__ == "__main__":
    sys.exit(main())
