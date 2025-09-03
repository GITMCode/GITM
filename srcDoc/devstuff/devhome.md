
# Developers

Here is some information which may be useful to developers.

## Documentation Notes

Some notes on the doc site:

- See the [markdown style guide](markdown_ref.md) for info on what features are
  available.
- The documentation is organized according to a table of contents listed in
  `GITM/mkdocs.yml`. So, to add new pages or organize things differently, look
  there.
- Within `mkdocs.yml`, paths are specified relative to `site_home` (cuttently
  set to `srcDoc/`)
- `index.md` is the "home page". The rest should be self-explanatory enough.

Writing docs can be tricky. Formatting is sometimes not what you would expect.
For example, we *need* four spaces to indent something, not just two. I tried
changing this & couldn't!

You can use Python to locally host the doc site which will auto-update when any
files are changed. One command to configure, one to host. When first starting,
either create a vitrual/conda environment for the docs, or don't, and from the
root of the repo run: `pip install -r srcDoc/requirements.txt`. To actually host
the docs, run (from `GITM/`):

```bash
mkdocs serve
```

And it should be auto-magic. open the link in your browser.

## Versioning

GITM now builds the code version into the executable automatically. This is sourced
from `src/.version` within `src/ModGITMVersion.f90`.

The version file is created when compiling, and checked every time the code is compiled.
If the number of files that are different from HEAD changes, the version is updated.
In practice this means that if the number of lines output from 
`git status --porcelain` changes, the version is updated.

The version has information on the last commit date, the branch, the last commit hash,
and the number of files changed since the last commit. This should provide enough
information to track down any kinds of bugs or to ensure that results can be traced
back to a certain state of the code.

The file is created in the last step of compiling by the shell script 
`share/Scripts/Makeversion.sh`.

If the `.git` folder does not exist (so the code was downloaded as a .zip from GitHub),
the version info is taken from `version.def` at the root folder. This should be updated
whenever the code is released.
