# Contributing Guidelines

Firstly, thank you for your interest in contributing to The Global
Ionosphere-Thermosphere Model! GITM is an open-source code and we welcome contributions
and feedback from users.

> This document serves to outline the requirements and best practices for those who
wish to contribute to GITM. Following these guidelines will ensure a positive
experience for all involved.

- [Contributing Guidelines](#contributing-guidelines)
  - [Opening an Issue](#opening-an-issue)
    - [Bug Report](#bug-report)
    - [Feature Requests](#feature-requests)
  - [Pull Requests](#pull-requests)
    - [Committing Changes](#committing-changes)
      - [Commit Styling](#commit-styling)
  - [Testing GITM](#testing-gitm)
  - [Code Formatting](#code-formatting)
    - [Automatic/Validating Code Formatting](#automaticvalidating-code-formatting)
    - [Style Guidelines \& `.fprettify.rc`](#style-guidelines--fprettifyrc)

---

## Opening an Issue

If you notice a bug, have a feature you would like to see, or have a problem with GITM,
the best practice is to [create an issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue#creating-an-issue-from-a-repository).

Within the Issues tab, there are two templates. Select which best describes the issue you
are creating and then modify the contents of the template to describe your issue.

### Bug Report

If you notice a problem with GITM and are unsure of how to fix it yourself, please
create a bug report detailing the problem. The template has several optional sections,
so read through them all before adding information to the wrong section.

- **Do not open an issue for problems installing GITM on a specific system.** Since
there is no way to predict all possible systems which GITM will be used on, details on
specific systems cannot be provided. Try the steps in the [README](../README.md) first,
then reach out to the development team if you need more help.
- Please include enough information to allow a maintainer to reproduce your bug,
only seeing the information in the bug report. Forgetting to attach input files or
only saying "the code did not run on system X", for example, will result in a delay.
- Check if another issue is open (or closed) with the same problem. Duplicate issues
may be deleted. To help prioritize issues, you can react or comment on an open issuue
to note that you are having the same problem.

### Feature Requests

Feature requests are most welcome! Please follow the template and provide as much detail
as you possibly can.

- Even if you are working on a feature yourself, creating an issue is not a bad idea!
Creating the issue will let other users know you are working on it and they may be able
to help. Please note in the issue that you are working on it and if you would like
help. Providing a link to the branch you are working on will help people find it easier
too.
- Not all features can be implemented. Sorry. We welcome as many ideas as possible, but
some things will be impossible, out-of-scope, or too time-consuming to implement.

## Pull Requests

 Before [forking the repo](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)
 and [creating a pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests),
 it is usually best to first open an issue to discuss the changes, or discuss your
 intended approach for solving the problem in the comments for an existing issue.
 But small changes may not warrant a new issue.

- **Pull Requests MUST go into the `develop` branch.** The `main` branch of GITM is meant
to remain stable and is only updated by releases, so new features and bug fixes should be
put into `develop` for testing
before going to `main`.
- **Submit pull requests for single changes at a time**. Large pull requests will both
take longer to review and may cause unintended problems for contributors working in
parallel. Please do not introduce changes to unrelated code within the same pull request
(example: refactoring one file and fixing a typo in a another file). If you
are working on a feature and notice an unrelated bug or typo, create a new branch from
`develop`, fix it, and submit the pull request separately. There is no harm in
submitting many small pull requests!
- **Document your changes**. If you are introducing a new feature, please include
a description of it and its uses in both the documentation (`srcDoc/` folder) and
in comments within the code.
- **Include new tests**. If you are adding a new feature, please add a sample `UAM.in`
file with this option enabled in the `srcTests/auto_test/` folder. We do not want future updates
breaking your hard work.
- **Follow existing code style guidelines**. More details can be found within the
[code_style](../code_style/) folder. We want to ensure the code remains readable and
consistent, so code formatting is tested upon the creation of a pull request. Incorrectly
formatted code will not be accepted.

### Committing Changes

1. First, fork the GITM repository. This will give you a clean slate to begin making changes.
2. Branch from `develop`. This ensures your changes are up-to-date with the latest patches of GITM.
3. Name your branch something that is descriptive!
4. Commit changes as they are completed, include a description for changes that are complicated or could break things.
5. Verify your formatting!

In short, make new branches for features `git checkout -b my_feature` and commit
often and push a little less often. Try to format then merge back to develop as soon as
you have something that works.

#### Commit Styling

The first line of the commit must be *at most* ~50 characters long and
should start with either:

- `FEAT:` For new feature.
- `BUG:` For bug fix.
- `MERGE:` For merging.
- `DOC:` For documentation update.
- `TEST:` For the addition or modification of tests.
- `STY:` For a style update (e.g., linting).
- `DEPREC:` Deprecate something, or removee a deprecated object.
- `REVERT:` Revert an earlier commit.
- `MAINT:` For maintenance such as refactoring, typos, etc.

The commit first line must be in *present* tense so that anyone
picking a commit hash can easily read what they are enabling. For more
information check out [conventional commit
messages](https://www.conventionalcommits.org/en/v1.0.0/).

For example,

*do:*

```gh
FEAT: Hydrostatic density implementation.
```

*don't:*

```gh
Implemented hydrostatic density. (feature)
```

## Testing GITM

GITM has a number of tests that are maintained and run automatically on every
release & pull-request. To run these yourself, run the script `run_all_tests.sh`
from within `srcTests/auto_test`. New tests can be added by simply creating
another UAM.in file, and will be run automatically if the file matches the
pattern `UAM.in.*.test`.

It is best practice to create tests as bugs are fixed. For example, if running
GITM in a certain configuration causes a crash, it is recommended to create a
test with this configuration which will help ensure the bug does not sneak back
in with future development.

## Code Formatting

GITM now uses a custom implementation of the
[`fprettify`](https://github.com/GITMCode/fprettify) library to format code. Pull
requests must be formatted prior to submission.

To install this version, (optionally create a new python environment and activate it, then) go the location you would like to download to. Then install fprettify with the commands:

```sh
git clone git@github.com:GITMCode/fprettify.git
cd fprettify
pip install .
```

Done! Run `fprettify -h` to see the options available from the command line.

More details on GITM's implementation of `fprettify` can be found on the [GITM/fprettify](https://github.com/GITMCode/fprettify) repository.

### Automatic/Validating Code Formatting

A python file has been made to aid in the formatting of code. To automatically format the entire GITM repository in-place, go to the root of the GITM repository and run

```sh
python srcPython/format_GITM.py -a
```

Upon creation of a pull request, the code will be checked for recommended changes. Before
submission, it may be a good idea to test it yourself with:

```sh
python srcPython/format_GITM.py -f
```

### Style Guidelines & `.fprettify.rc`

The code guidelines are listed in `.fprettify.rc`, which is automatically sourced from by
`code_formatter.py`. This file is included within GITMCode/GITM and will be found
automatically when running from the root of the GITM repository. If you call `format_GITM.py`
from another directory, `src/` for example, there is a chance that configuration file will not
be found. The path to your `.fprettify.rc` file can be specified with the `--fpretty_config`
flag when running `format_GITM.py`, or you can copy this file in your `$HOME` directory, where
it will always be found.

The config file specifies the following requirements:

1. Maximum line length of 120 characters. 88 is recommended but not always possible.
2. Indent with 2 spaces. No tabs!
3. Check all files, except: anything compiled, src/Venus.f90, old fortran files, makefiles, MPI Templates
4. A single space after:  operators (`for`, `if`, etc.), plus & minus
5. No spaces after: intrinsics (`write(`, *not* `write (` for example), multiple & divide.

Some recommended guidelines are:

1. All lines shorter than 88 characters. Use a `&` to break long lines
2. Subroutines shorter than 250 lines. Refactoring long subroutines can help a lot with readability.
