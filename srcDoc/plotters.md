
# Plotters

## gitm_plot_simple.py

This code will make simple rectangular plots in:

- altitude (longitude vs latitude): '-cut=alt -alt=value'
- longitude (latitude vs altitude): '-cut=lon -lon=value'
- latitude (longitude vs altitude): '-cut=lat -lat=value'

You can set the variable to plot as a number. To see the variable mapping between numbers and names, you can run the code with '-list'. You can set the max and min value for the colorbar if you would like, but setting -maxi and -mini.

    usage: gitm_plot_simple.py [-h] [-list] [-var VAR] [-cut cut] [-alt alt] [-lat lat] [-lon lon]
                               [-log] [-mini MINI] [-maxi MAXI]
                               filelist [filelist ...]

    Plot Aether / GITM model results

    positional arguments:
        filelist    list files to use for generating plots

    options:
        -h, --help  show this help message and exit
        -list       list variables in file
        -var VAR    variable to plot (number)
        -cut cut    alt,lat,lon : which cut you would like
        -alt alt    altitude : alt in km (closest)
        -lat lat    latitude : latitude in degrees (closest)
        -lon lon    longitude in degrees (closest)
        -log        plot the log of the variable
        -mini MINI  manually set the minimum value for the plots
        -maxi MAXI  manually set the maxiumum value for the plots

## gitm_plot.py

This code is a bit more sophisticated, adding a bunch of additional functionality:

- It can make cuts in alt, lat, and lon, as above.
- The alt cut plots add circles for the northern and southern polar regions.  These can be removed with the flag '-nopole'. Additionally, if you only want a polar plot, you can use the flag '-north' or '-south'.
- TEC and O/N2 can be plotted by using the flags '-tec' or '-on2'.  You can't set the both at the same time.
- You can overplot winds by setting the flag '-winds'. The code tries to figure out whether you want neutral winds or ion drifts by looking at the variable number and seeing if the main variable to plot is a neutral or ion state. You can set the spacing between the wind vectors by setting '-nstep=value'.
- You can subtract a series of files from another series of files.  If the '-backdir=/path/to/baseline/run' is set, it will look into this directory and try to find matching files to subtract. If you set '-percent', it computes the percentage difference.
- If a series of images is made (by feeding a bunch of files), you can make a movie by setting the flag '-movie'. You can set the framerate ('-rate') and type of movie ('-mkv', '-mp4', '-gif'). You need ffmpeg to actually use this feature.
- The code can integrate or calculate the average of the slice(s) and then plot those out to a file by setting '-timeplot'. Additionally, a logfile can be created by setting '-timefile=log_file_name' to output those means or integrals out to the log file.

    usage: thermo_plot.py [-h] [-list] [-timeplot] [-mean] [-timefile TIMEFILE] [-label]
                          [-color {default,red}] [-var VAR] [-cut cut] [-ext {png,jpg,pdf}]
                          [-hiwind HIWIND] [-winds] [-nstep nstep] [-nopole] [-north] [-south]
                          [-alt alt] [-lat lat] [-lon lon] [-alog] [-IsLog] [-mini MINI] [-maxi MAXI]
                          [-percent] [-diff] [-backdir backdir] [-mkv] [-mp4] [-gif] [-movie] [-tec]
                          [-on2] [-rate RATE]
                        filelist [filelist ...]

    Plot Aether / GITM model results

    positional arguments:
        filelist              list files to use for generating plots

    options:
        -h, --help            show this help message and exit
        -list                 list variables in file
        -timeplot             Plot integrated (or mean) value vs. time
        -mean                 Plot mean value instead of integrated value
        -timefile TIMEFILE    output filename for timeline file
        -label                Add label (e.g., (a), (b)..) to title
        -color {default,red}  set color bar
        -var VAR              variable to plot (number)
        -cut cut              alt,lat,lon : which cut you would like
        -ext {png,jpg,pdf}    plot type file extention
        -hiwind HIWIND        HIWIND file to plot location
        -winds                overplot winds
        -nstep nstep          number of steps between wind quivers
        -nopole               dont plot polar regions
        -north                only plot northern hemisphere results
        -south                only plot southern hemisphere results
        -alt alt              altitude : alt in km (closest)
        -lat lat              latitude : latitude in degrees (closest)
        -lon lon              longitude in degrees (closest)
        -alog                 plot the log of the variable
        -IsLog                plot the log of the variable
        -mini MINI            manually set the minimum value for the plots
        -maxi MAXI            manually set the maxiumum value for the plots
        -percent              plot percentage difference of files
        -diff                 plot difference of files (2 files needed)
        -backdir backdir      Subtract files in this directory
        -mkv                  movie format = mkv
        -mp4                  movie format = mp4
        -gif                  movie format = gif
        -movie                Make a movie out of results
        -tec                  plot total electron content (TEC)
        -on2                  plot O/N2 ratio
        -rate RATE            framerate for movie

In order to run this code, you need to install aetherpy, which can be obtained [here](https://github.com/aaronjridley/aetherpy). The 'develop' branch is fully supported, so please switch to that branch ('git checkout develop') and follow the README.md file to install it.  It should be straightforward to install.

## Aetherpy

Aetherpy is a python package for plotting Aether model results, but we have modified it so that it will support GITM files also. Aetherpy can be obtained [here](https://github.com/aaronjridley/aetherpy). The 'develop' branch is fully supported, so please switch to that branch ('git checkout develop') and follow the README.md file to install it.  It should be straightforward to install.

There are plotters in aetherpy which are useful but somewhat coomplicated. We are working on these to make them less complicated and more functional.

## PyITM

We are working on a new repository for plotting IT modeling results, and have put it [here](https://github.com/GITMCode/PyITM).  This is relatively new.

This code base is aimed at getting people ramped up to plotting and analysis quickly, so we are aiming for the codes to be relatively simple and easy to use. In addition, we are trying to make functions short and reusable so that they can be taken from these codes and used again (for example, the 'gitm_plot_simple.py' code is almost completely stolen from these functions and was created in about an hour from this framework.)

## Data-Model Comparison Tools

Some codes that have been provided to allow data-model comparisons to be made include the following below.  For each, to get more information, run them with a '-h'.  We don't provide this data to you - it is up to you to download and put this data somewhere that the codes can find it.  For almost all of these files, the path to the data is hardcoded into the file, so you will have to edit the file and replace the location of the data to whereever you install the data.

- thermo_goce.py: This code will compare GITM 3DALL files to GOCE measurements of densities and winds. Since GOCE was only available between late 2009 - 2013, if it doesn't find GOCE files, it will search for CHAMP files.  If those are not available, it will search for SWARM files. To download (for example) CHAMP data, you can use the command:
```bash
wget --no-parent -r --cut-dirs=3 ftp://thermosphere.tudelft.nl/version_01/CHAMP_data/
```
Other data types like this can be found on the same ftp site.

- thermo_guvi.py: This code will compare the GITM 3DALL files to (post 2006) GUVI O/N2. 

- thermo_vista.py: This code will read in and do some crude comparisons between GITM and VISTA TEC.  It is a work in progress.

- thermo_hiwind.py: This code will compare the HIWIND FPI results (for extremely specific times in 2011 and 2018) to 3DALL results.

