#!/usr/bin/env python

import numpy as np
import argparse, sys
import datetime
from itertools import zip_longest

## DEFAULTS
RELTOL = 1e-6
ATOL = 1e-3
START_STRING = "#START"


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


def assess_diffs(name, arr1, arr2, times=None, reltol=RELTOL, atol=ATOL):
    """
    returns a string with info about how/where two arrays differ
    """

    # Boolean array of compliant values
    equal_ma = np.isclose(arr1, arr2, rtol=reltol, atol=atol, equal_nan=False)
    nDifferent = len(equal_ma) - np.count_nonzero(equal_ma)
    firstDifferent = int(np.argmin(equal_ma))

    msg = f">>> {name}: ({nDifferent} / {len(equal_ma)}) values differ, first at "

    if times is not None:
        msg += f"+{(times[firstDifferent] - times[0]).total_seconds():g} s, or "
    msg += f"index {firstDifferent}"

    return msg


def main(file1, file2, reltol=RELTOL, atol=ATOL, verbose=False):

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
        s += f"   {'file 1':<{w}}file 2"
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
                f"\nfile 1 has {nLines1} lines, file 2 has {nLines2}")
        # Only compare over the overlap
        minLines = min(nLines1, nLines2)
        for log in (log1, log2):
            for k in log:
                log[k] = log[k][:minLines]

    # Then check the numerics, skipping any cols necessary
    skip_cols = ["time"]
    failed = []
    for name in log1:
        if name in skip_cols or name not in log2:
            continue
        if np.allclose(
            log1[name], log2[name], rtol=reltol, atol=atol, equal_nan=False
        ):
            if verbose:
                print(f"   {name}: ok")
        else:
            s = assess_diffs(
                name, log1[name], log2[name],
                times=log1.get("time"), reltol=reltol, atol=atol)
            failed.append(s)
    errors.extend(failed)
    if not errors:
        sys.exit(0)

    # Name the files once, here, so the lines above stay short
    print("\n" + "=" * 70 + f"\n   file 1:  {file1}\n   file 2:  {file2}\n" + "-" * 70)
    print("Test did not pass!")
    for s in errors:
        print("\n" + s)
    sys.exit(1)


if __name__ == "__main__":
    args = parse_args()
    main(
        args.file1,
        args.file2,
        reltol=args.rtol,
        atol=args.atol,
        verbose=args.verbose,
    )
