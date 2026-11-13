#!/usr/bin/env python3
"""Generate the per-type variable tables in srcDoc/outputs.md.

GITM writes UA/data/output_schema.csv at the first output of every run: one
row per variable of every output type (type, variable, units, long_name),
whatever the backend and whether or not the type was requested. A copy is
kept in srcDoc/outputs/output_schema.csv, and the tables in srcDoc/outputs.md
are generated from that copy.

After changing an output type, run any GITM case and then:

    python srcPython/generate_output_docs.py <run>/UA/data/output_schema.csv

which replaces the copy in srcDoc/outputs/ and rewrites the tables between
the BEGIN/END markers in srcDoc/outputs.md. With no argument the tables are
rebuilt from the existing copy.

    python srcPython/generate_output_docs.py --check [<run>/UA/data/output_schema.csv]

changes nothing and exits 1 if the docs are stale: the given file differs
from the copy, the tables differ from what the copy generates, or a type is
missing from the hand-written part of outputs.md. run_all_tests.sh runs this
after the tests. Needs only the standard library.
"""

import argparse
import csv
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
SCHEMA = os.path.join(REPO, "srcDoc", "outputs", "output_schema.csv")
DOC = os.path.join(REPO, "srcDoc", "outputs.md")

BEGIN = "<!-- BEGIN GENERATED TABLES: edit with srcPython/generate_output_docs.py, not by hand -->"
END = "<!-- END GENERATED TABLES -->"

# Valid type codes that stop the run when requested (their fills are not
# wired to the GLOW emissions), so they get no table. outputs.md says why.
NO_TABLE = {"3DGLO", "1DGLO"}


def read_rows(path):
    """The schema rows of a csv file, skipping its '#' comment lines."""
    with open(path, newline="") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines))


def types_in(rows):
    """Type codes in table order: 3D first, then 2D, 1D, 0D; A-Z within each."""
    return sorted({r["type"] for r in rows}, key=lambda t: (-int(t[0]), t))


def cell(text):
    """A string made safe for a markdown table cell."""
    return text.strip().replace("|", "\\|")


def make_tables(rows):
    out = []
    for t in types_in(rows):
        if t in NO_TABLE:
            continue
        out += [f"### {t}", "", "| Variable | Units | Description |", "|---|---|---|"]
        out += [f"| {cell(r['variable'])} | {cell(r['units'])} | {cell(r['long_name'])} |"
                for r in rows if r["type"] == t]
        out.append("")
    return "\n".join(out)


def split_doc(text):
    """outputs.md as (before, generated, after), split at the two markers."""
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        sys.exit(f"{DOC}: expected exactly one BEGIN and one END marker line")
    before, rest = text.split(BEGIN)
    generated, after = rest.split(END)
    return before, generated, after


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("schema", nargs="?",
                        help="an output_schema.csv written by a GITM run")
    parser.add_argument("--check", action="store_true",
                        help="report stale docs and exit 1; change nothing")
    args = parser.parse_args()

    rows = read_rows(args.schema or SCHEMA)
    if not rows:
        sys.exit(f"no schema rows in {args.schema or SCHEMA}")
    with open(DOC) as f:
        before, generated, after = split_doc(f.read())
    tables = "\n\n" + make_tables(rows) + "\n"

    # A type the hand-written text never names has no description anywhere.
    undescribed = [t for t in types_in(rows) if f"`{t}`" not in before + after]

    if args.check:
        problems = []
        if args.schema and read_rows(SCHEMA) != rows:
            problems.append(f"{os.path.relpath(SCHEMA, REPO)} does not match {args.schema}")
        if generated != tables:
            problems.append("the variable tables in srcDoc/outputs.md are out of date")
        if undescribed:
            problems.append("srcDoc/outputs.md does not describe: " + ", ".join(undescribed))
        if problems:
            print("Output docs are stale:")
            for p in problems:
                print("  - " + p)
            print("To update, from the top of the repository:")
            source = os.path.abspath(args.schema) if args.schema else ""
            print(f"  python srcPython/generate_output_docs.py {source}".rstrip())
            sys.exit(1)
        print(f"Output docs match the schema ({len(types_in(rows))} types).")
        return

    if args.schema and os.path.realpath(args.schema) != SCHEMA:
        with open(args.schema) as src, open(SCHEMA, "w") as dst:
            dst.write(src.read())
        print("updated", os.path.relpath(SCHEMA, REPO))
    with open(DOC, "w") as f:
        f.write(before + BEGIN + tables + END + after)
    print("updated", os.path.relpath(DOC, REPO))
    if undescribed:
        print("Not described in srcDoc/outputs.md yet (add them to the type list): "
              + ", ".join(undescribed))


if __name__ == "__main__":
    main()
