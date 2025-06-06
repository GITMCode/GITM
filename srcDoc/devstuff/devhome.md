
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

## To Do list for documentation

Feel free to add to this

- [ ] Section on FAQ?
- [ ] More on internals:
    - [ ] Electrodynamics section!
    - [ ] Better info in common inputs?
- [ ] GitHub walkthrough in dev section
- [ ] some way to auto-update outline?
- [ ] more more more!
