# Introduction

Global Ionosphere Thermosphere Model

## Welcome

Hi. Thanks for visiting this website. Here is some info on how to use GITM.

Please see the [quick start](quick_start.md) page for more info on how to
download, configure, and install GITM.


---

## Notes for Developers

Some notes on the doc site:

- See the [markdown style guide](markdown_ref.md) for info on what features are
  available.
- The documentation, presently, looks in `srcDoc/` for any markdown files and
  adds them to this site alphabetically. This can be changed later, but filling
  the docs out initially will be easier if pages are automatically added. But
  that's why there is a page for the development team.
- `index.md` is the "home page". The rest should be self-explanatory

Writing docs can be tricky. Formatting is sometimes not what you would expect.
For example, we *need* four spaces to indent something, not just two. I tried
changing this & couldn't!

You can use Python to locally host the doc site which will auto-update when any
files are changed. One command to configure, one to host. When first starting,
either create a vitrual/conda environment for the docs, or don't, and from the
root of the repo run: `pip install -r srcDoc/requirements.txt`. To actually host
the docs, run:

```bash
mkdocs serve
```

And it should be auto-magic. open the link in your browser.

!!! warning
    Make sure to run this from the root of the repository! Otherwise it won't
    work!

### To Do list for documentation

Feel free to add to this

- [x] Move over the existing latex manual pages
    - [ ] Verify that info is correct, relevant, and necessary (ALB in progress)
    - [ ] Organize that information (ALB in progress)
    - [ ] Fix links from latex manual pages
- [ ] Keep filling out information!


### Draft outline for docs

This is off the cuff & is open to feedback. Just change it to what you want. A
few sections here are repetitive and don't need to be.

- welcome page
    - Quick overview of what GITM is
    - installation & config (super basic)
- Quick start
    - more detailed config options
    - how to actually run the code
- Grid
    - cells vs blocks
    - how to get to the resolution you want
- Inputs - common
    - Important and/or frequently changed UAM options
    - more details on files necessary for some options
- Inputs - all
    - long page with all possible options that GITM can check
- Outputs
    - each output type's variables?
  - FAQ
    - common issues:
        - indices times are wrong
        - running with a different compiler than compiled with?
        - how do I configure for X system? (known good configs)
    - where can I find example input files? (GITMCode/GITM_Input_files)
    - where can I find example outputs (CCMC?)
    - Can you add this feature? no. You can, or fill out a feature request.
- Electrodynamics
- Chemistry
