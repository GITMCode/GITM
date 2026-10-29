#!/usr/bin/env python

import numpy as np
import argparse, sys
import datetime
import difflib
from itertools import zip_longest

## DEFAULTS
RELTOL = 1e-6
ATOL = 1e-3
TIMETOL = 0.25  # seconds
START_STRING = "#START"
NROWS = 5  # rows printed at the head and at the tail


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare two GITM logfiles column-by-column.",
        epilog="Exits 0 if the logs match, 1 otherwise.",
    )
    parser.add_argument("file1", help="reference logfile")
    parser.add_argument("file2", help="logfile to test")
    parser.add_argument(
        "-r",
        "--rtol",
        type=float,
        default=RELTOL,
        help=f"relative tolerance (default {RELTOL:g})",
    )
    parser.add_argument(
        "-a",
        "--atol",
        type=float,
        default=ATOL,
        help=f"absolute tolerance (default {ATOL:g})",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="report every column, not just failures",
    )
    parser.add_argument(
        "-n",
        "--rows",
        type=int,
        default=NROWS,
        help=f"rows to print at head and tail, 0 for all (default {NROWS})",
    )
    # Accepted for compatibility with DiffNum.pl call sites; no text channel here.
    parser.add_argument(
        "-b", "--blank", "--ignore-blanks", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-t", "--text", "--ignore-text", action="store_true", help=argparse.SUPPRESS
    )

    return parser.parse_args()


def read_file(fname, headerText=START_STRING, parse_datetimes=True):
    """A simple GITM logfile reader

    Inputs
    ------
      fname (str): path to logfile
      headerText (str): what appears in the line before data starts
      parse_datetimes (bool): Transform (year, month, etc.) to a python datetime?
                              default=True

    Returns
    -------
      (dict): with keys of the column (or time) & values as np arrays

    """

    outData = {}
    didStart = False
    isFirstLine = False

    with open(fname, "r") as file:
        allLines = file.readlines()
        for iLine, line in enumerate(allLines):
            if headerText in line:
                didStart = True
                isFirstLine = True
                continue

            if isFirstLine:
                for word in line.strip().split():
                    outData[word] = []
                isFirstLine = False
                continue

            if didStart:
                vals = line.split()
                if not vals:
                    continue
                if len(vals) != len(outData):
                    # a killed run leaves a partial final line
                    print(f">> {fname}: skipping line {iLine + 1}, "
                          f"got {len(vals)} of {len(outData)} columns")
                    continue
                for key, val in zip(outData.keys(), vals):
                    outData[key].append(val)

        if not didStart:
            raise ValueError(
                f"Header '{headerText}' not found in logfile! Could not parse"
            )

    if parse_datetimes:
        # Try the default datetime keys:
        dt_vars = ["yyyy", "mm", "dd", "HH", "MM", "SS", "ms"]
        if not all([v in outData.keys() for v in dt_vars]):
            print(
                ">> Default time keys not found by DiffNum.py. Did the logfile format change?"
            )
            # Try the [1-8] variables in the header
            # should be "(iStep, yyyy, mm, dd, HH, MM, SS, MS, [...])"
            dt_vars = list(outData.keys())[1:8]

        # parse & put into a numpy array
        times = np.array([
            datetime.datetime(*[int(outData[k][i]) for k in dt_vars[:6]])
            + datetime.timedelta(milliseconds=int(outData[dt_vars[6]][i]))
            for i in range(len(outData[dt_vars[0]]))
        ])

        for k in dt_vars:
            outData.pop(k)
        outData["time"] = times
    for k in outData:
        if k != "iStep" and k != "time":
            outData[k] = np.array(outData[k], dtype=float)
        elif k == "iStep":
            outData[k] = np.array(outData[k], dtype=int)

    return outData


def read_header(fname, headerText=START_STRING):
    """Return the lines above the data block (the echoed inputs)"""
    header = []
    with open(fname, "r") as file:
        for line in file:
            if headerText in line:
                break
            header.append(line.rstrip())
    return header


def pick_rows(first, nTotal, n):
    """Head (from the first difference, or the top if it is near) and tail"""
    if n <= 0:
        return list(range(nTotal))
    top = 0 if first < n else first
    head = list(range(top, min(top + n, nTotal)))
    tail = list(range(max(nTotal - n, head[-1] + 1), nTotal))
    return head + tail


def value_table(names, log1, log2, bad, rows, times=None):
    """One line per (column, file), one field per printed row.

    Time runs across the page, so any number of failing columns fits;
    a gap in the printed rows shows up as a '...' field.
    """
    def num(v):
        # the logs write %.3f; keep that, unless it would hide the value
        if v == 0 or 1e-3 <= abs(v) < 1e6:
            return f"{v:.3f}"
        return f"{v:.4e}"

    # Field list: a row index, or None for a gap
    fields = []
    for i in rows:
        if fields and i != fields[-1] + 1:
            fields.append(None)
        fields.append(i)

    steps = log1.get("iStep", np.arange(len(log1[names[0]])))
    cells = {(k, f): [num(log[k][i]) if i is not None else "..." for i in fields]
             for k in names for f, log in (("ref", log1), ("test", log2))}
    tops = [str(steps[i]) if i is not None else "..." for i in fields]
    hms = [times[i].strftime("%H:%M:%S") if times is not None and i is not None
           else "" for i in fields]
    w = max(len(c) for v in list(cells.values()) + [tops, hms] for c in v) + 1
    lw = max(len(k) for k in names + ["iStep"]) + 1

    def line(label, tag, vals):
        return f"{label:<{lw}}{tag:<5}" + "".join(f"{v:>{w}}" for v in vals)

    out = [line("iStep", "", tops)]
    if times is not None:
        out.append(line("", "", hms))
    for k in names:
        out.append("")
        out.append(line(k, "ref", cells[(k, "ref")]))
        out.append(line(f"({np.count_nonzero(bad[k])})", "test", cells[(k, "test")]))
    return "\n".join(out)


def main(file1, file2, reltol=RELTOL, atol=ATOL, verbose=False, nRows=NROWS):

    # Run every check, then fail at the end if any did not pass
    errors = []

    # Read logs, parsing datetimes
    log1 = read_file(file1)
    log2 = read_file(file2)

    # First make sure the column names match. Compared by name below, so
    # a reordering is harmless -- only unmatched names are worth reporting.
    only1 = [k for k in log1 if k not in log2]
    only2 = [k for k in log2 if k not in log1]
    if only1 or only2:
        w = max(len(k) for k in only1 + only2) + 3
        s = "Unmatched column names, not compared:\n"
        s += f"   {'ref':<{w}}test"
        for a, b in zip_longest(only1, only2, fillvalue="-"):
            s += f"\n   {a:<{w}}{b}"
        errors.append(s)

    # Then make sure the logfiles are the same length
    nLines1 = len(next(iter(log1.values())))
    nLines2 = len(next(iter(log2.values())))
    if nLines1 != nLines2:
        if abs(nLines1 - nLines2) == 1:
            # Sometimes the logfiles only differ in the last line. 
            # Do not cause the test to fail, but warn that things are different
            print("Number of lines in the two files are off by 1")
        else:
            errors.append(
                f"\nref has {nLines1} lines, test has {nLines2}")
        # Only compare over the overlap
        minLines = min(nLines1, nLines2)
        for log in (log1, log2):
            for k in log:
                log[k] = log[k][:minLines]

    # Then make sure the runs are in lockstep. Comparing columns row by row
    # only means something if both rows are the same moment of the simulation.
    t1, t2 = log1.get("time"), log2.get("time")
    if t1 is not None and t2 is not None:
        n = min(len(t1), len(t2))
        drift = np.array([(t2[i] - t1[i]).total_seconds() for i in range(n)])
        bad = np.abs(drift) > TIMETOL
        if bad.any():
            i = int(np.argmax(bad))
            errors.append(
                f"\nTimes differ: {np.count_nonzero(bad)} / {n} rows "
                f"\n   first at row {i}:   ref {t1[i]}   test {t2[i]}")

    # Then check the numerics, skipping any cols necessary
    skip_cols = ["time"]
    compared = [k for k in log1 if k not in skip_cols and k in log2]
    bad = {k: ~np.isclose(log1[k], log2[k], rtol=reltol, atol=atol,
                          equal_nan=False) for k in compared}
    failed = [k for k in compared if bad[k].any()]
    if verbose:
        for k in compared:
            if k not in failed:
                print(f"   {k}: ok")
    if not errors and not failed:
        sys.exit(0)

    # Name the files once, here, so the lines above stay short
    print("\n" + "=" * 70 + f"\n   ref:   {file1}\n   test:  {file2}\n" + "-" * 70)
    print(f"Test did not pass!   (rtol {reltol:g}, atol {atol:g})")
    for s in errors:
        print(s)

    if failed:
        times = log1.get("time")
        nTotal = len(log1[failed[0]])
        anyBad = np.any([bad[k] for k in failed], axis=0)
        first = int(np.argmax(anyBad))
        where = f"iStep {log1['iStep'][first]}" if "iStep" in log1 else f"row {first}"
        if times is not None:
            where += f" (+{(times[first] - times[0]).total_seconds():g} s)"
        print(f"\n{len(failed)} of {len(compared)} columns differ, in "
              f"{np.count_nonzero(anyBad)} of {nTotal} rows; first at {where}."
              f"\n(n) under a name = rows that column differs in."
              f"  -n sets rows shown at head and tail, -n 0 all.\n")
        print(value_table(failed, log1, log2, bad,
                          pick_rows(first, nTotal, nRows), times) + "\n")

    # The echoed inputs often say why: report them, but never fail on them
    h1, h2 = read_header(file1), read_header(file2)
    if h1 != h2:
        diff = list(difflib.unified_diff(h1, h2, "ref", "test",
                                         n=0, lineterm=""))
        print("Header (echoed inputs) also differs:")
        for line in diff[2:]:
            if line[0] in "-+":
                print(f"   {'ref only: ' if line[0] == '-' else 'test only:'} {line[1:]}")
    sys.exit(1)


if __name__ == "__main__":
    args = parse_args()
    main(
        args.file1,
        args.file2,
        reltol=args.rtol,
        atol=args.atol,
        verbose=args.verbose,
        nRows=args.rows,
    )
