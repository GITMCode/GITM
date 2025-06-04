#!/usr/bin/env python3

"""
I can't write perl, so this is a python wrapper around src/inputs_to_tex.pl


It requires pandoc to bt installed since it converts the .tex file to .md

Does not need to be re-ran unless set_inputs.f90 is changed.


usage:

./set_inputs.py

or 

python set_inputs.py

Must be within srcDoc!

"""

import os, subprocess

# run src/inputs_to_tex.pl
os.chdir('../src/')
tex2pl = subprocess.run("./inputs_to_tex.pl", shell=True, capture_output=True)

if tex2pl.returncode != 0:
    print(tex2pl.stderr.decode())
    print(tex2pl.stdout.decode())
    raise ValueError("Could not run src/inputs_to_tex.pl")

pandocing = subprocess.run("pandoc set_inputs.tex -o ../srcDoc/set_inputs.md",
                            shell=True, capture_output=True)

os.remove('set_inputs.tex')

if pandocing.returncode != 0:
    print(pandocing)
    raise ValueError("Could not run pandic. Is it installed??")


# Go back to srcdoc
# Rearrange set_inputs.md to have markdown headers
os.chdir("../srcDoc")

outlines = ['# All Inputs\n\n',
            '<!-- This file is automatically made by set_inputs.py -->\n\n'
            ]
this_sec = ['']
started_indent = False


with open('set_inputs.md') as f:
    alllines = f.readlines()

    for aline in alllines:
        this_sec.append(aline)
        if aline[:5] == '    #':
            this_sec[0] = aline.replace('    #', '## ') + '\n'
            started_indent = True
        if aline[:4] != '    ' and started_indent:
            outlines.extend(this_sec)
            this_sec = ['']
            started_indent = False

os.remove('set_inputs.md', )
with open('set_inputs.md', 'x') as f:
    for i in outlines:
        f.write(i)


print("all done")
