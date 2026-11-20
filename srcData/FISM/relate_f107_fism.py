#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from datetime import datetime
import numpy as np
import os, glob
from pyitm.fileio import fism
from pyitm.plotting import line_plots, axes

iLowNorm = 15
iLowCutoff = 25
iHighNorm = 85
iHighCutoff = 75
nTrials = 1000
nWaves = 59

def read_waves(euvCsv):
    nSkip = 5
    with open(euvCsv, 'r') as fpin:
        allLines = fpin.readlines()
        values = allLines[0].strip().split(',')
        shorts = np.flip(np.array(values[nSkip:-1]).astype(float))/10.0
        values = allLines[1].strip().split(',')
        longs = np.flip(np.array(values[nSkip:-1]).astype(float))/10.0
        fpin.close()
    return shorts, longs

#-----------------------------------------------------------------------------
# Custom error calculator
#-----------------------------------------------------------------------------

def calc_error(quantity, model, donorm = True):
    # this will be a bit complicated, since we want to capture highs and lows:
    # lows first:
    if (donorm):
        norm = np.percentile(quantity, 10)
    else:
        norm = 1.0
    highest = np.percentile(quantity, 15)
    qLow = quantity[quantity < highest]
    mLow = model[quantity < highest]
    errorLow = np.mean(np.abs(qLow - mLow) / norm)
    
    # highs next:
    if (donorm):
        norm = np.percentile(quantity, 90)
    else:
        norm = 1.0
    lowest = np.percentile(quantity, 60)
    qHigh = quantity[quantity > lowest]
    mHigh = model[quantity > lowest]
    errorHigh = np.mean(np.abs(qHigh - mHigh) / norm)

    # overall:
    if (donorm):
        norm = np.median(quantity)
    else:
        norm = 1.0
    errorOverall = np.median(np.abs(quantity - model) / norm)

    errorTotal = 5.0 * errorLow + 7.5 * errorHigh + 2.0 * errorOverall
    
    return errorTotal

#-----------------------------------------------------------------------------
# use monte carlo method to find best fit parameters for the model
#-----------------------------------------------------------------------------

def fit_model(f107, f107a, ir, donorm = True):

    p107 = (f107 + f107a)/2.0

    rF107 = np.percentile(p107, 90) - np.percentile(p107, 10)
    rIr = np.percentile(ir, 90) - np.percentile(ir,10)

    slopeF = rIr / rF107 / 2
    slopeA = slopeF
    inter = np.percentile(ir,10) - np.percentile(p107, 10) * slopeF * 2

    # Make an initial guess, which is going to be pretty bad:
    slopeFsave = slopeF
    slopeAsave = slopeA
    intersave = inter
    if (intersave < 0):
        intersave = np.percentile(ir, 10)
    pFsave = 1.0
    pAsave = 1.0
    
    model = slopeFsave * (f107**pFsave) + slopeAsave * (f107a**pAsave) + intersave
    model[model < 0.0] = 0.0
    rms = calc_error(ir, model, donorm = donorm)
    print('Initial Error : ', rms, slopeFsave, pFsave)
    rmssave = rms
    for iT in range(5):
        nT = nTrials * (iT + 1)
        print('Loop : ', iT)
        facs = [0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.01, 0.005]
        for fac in facs:
            print(' -> fac : ', fac)
            slopeFs = np.random.normal(loc=slopeFsave, scale=np.abs(slopeFsave) * fac, size=nT)
            slopeAs = np.random.normal(loc=slopeAsave, scale=np.abs(slopeAsave) * fac, size=nT)
            inters = np.random.normal(loc=intersave, scale=np.abs(intersave) * fac, size=nT)

            pFs = np.random.normal(loc=pFsave, scale=np.abs(pFsave) * fac*0.25, size=nT)
            pAs = np.random.normal(loc=pAsave, scale=np.abs(pAsave) * fac*0.25, size=nT)

            for i in range(nT):
                model = \
                  slopeFs[i] * (f107**pFs[i]) + \
                  slopeAs[i] * (f107a**pAs[i]) + inters[i]
                model[model < 0.0] = 0.0
                rms = calc_error(ir, model, donorm = donorm)
                if (rms < rmssave):
                    print('Better solution found! ', rms, np.mean(model), intersave)
                    slopeFsave = slopeFs[i]
                    slopeAsave = slopeAs[i]
                    intersave = inters[i]
                    pFsave = pFs[i]
                    pAsave = pAs[i]
                    rmssave = rms
        
    return slopeFsave, slopeAsave, intersave, pFsave, pAsave

#-----------------------------------------------------------------------------
# use monte carlo method to find best fit parameters for the model
#-----------------------------------------------------------------------------

def fit_model2(f107, f107a, ir, donorm = True):

    p107 = (f107 + f107a)/2.0

    rF107 = np.percentile(p107, 90) - np.percentile(p107, 10)
    rIr = np.percentile(ir, 90) - np.percentile(ir,10)

    slopeF = rIr / rF107 / 2
    slopeA = slopeF
    inter = np.percentile(ir,10) - np.percentile(p107, 10) * slopeF * 2

    # Make an initial guess, which is going to be pretty bad:
    slopeFsave = slopeF
    slopeAsave = slopeA
    intersave = inter
    if (intersave < 0):
        intersave = np.percentile(ir, 10)
    pFsave = 1.0
    pAsave = 1.0
    
    model = slopeFsave * (f107 - f107a) + slopeAsave * (f107a**pAsave) + intersave
    model[model < 0.0] = 0.0
    rms = calc_error(ir, model, donorm = donorm)
    print('Initial Error : ', rms, slopeFsave, pFsave)
    rmssave = rms
    for iT in range(5):
        nT = nTrials * (iT + 1)
        print('Loop : ', iT)
        facs = [0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.01, 0.005]
        for fac in facs:
            print(' -> fac : ', fac)
            slopeFs = np.random.normal(loc=slopeFsave, scale=np.abs(slopeFsave) * fac, size=nT)
            slopeAs = np.random.normal(loc=slopeAsave, scale=np.abs(slopeAsave) * fac, size=nT)
            inters = np.random.normal(loc=intersave, scale=np.abs(intersave) * fac, size=nT)

            pAs = np.random.normal(loc=pAsave, scale=np.abs(pAsave) * fac*0.25, size=nT)

            for i in range(nT):
                model = \
                  slopeFs[i] * (f107 - f107a) + \
                  slopeAs[i] * (f107a**pAs[i]) + inters[i]
                model[model < 0.0] = 0.0
                rms = calc_error(ir, model, donorm = donorm)
                if (rms < rmssave):
                    print('Better solution found! ', rms, np.mean(model), intersave)
                    slopeFsave = slopeFs[i]
                    slopeAsave = slopeAs[i]
                    intersave = inters[i]
                    pAsave = pAs[i]
                    rmssave = rms
        
    return slopeFsave, slopeAsave, intersave, pFsave, pAsave


allFiles = sorted(glob.glob('fismflux_daily_*.dat'))
f107File = '../f107.txt'

fismData = fism.read_all_fism_files(allFiles)
f107Data = fism.read_f107_file(f107File)

print(fismData['times'], fismData['irradiances'][:,0])
print(np.shape(fismData['irradiances'][:,:]))

nFism = len(fismData['times'])
f107 = np.zeros(nFism)
f107a = np.zeros(nFism)
# assume that both fism and f107 are sorted
iF107 = 0
for iFism in range(nFism):
    tFism = fismData['times'][iFism]
    while ((f107Data['times'][iF107] - tFism).total_seconds() < 0):
        iF107 = iF107 + 1
    dtb = (f107Data['times'][iF107]-f107Data['times'][iF107-1]).total_seconds()
    dtt = (tFism - f107Data['times'][iF107-1]).total_seconds()
    x = dtt / dtb
    f107[iFism] = (1.0 - x) * f107Data['f107'][iF107-1] + \
      (x) * f107Data['f107'][iF107]
    f107a[iFism] = (1.0 - x) * f107Data['f107a'][iF107-1] + \
      (x) * f107Data['f107a'][iF107]



fileout1 = 'euv_model1.txt'
fileout2 = 'euv_model2.txt'
euvCsv = 'euv_59_v2.csv'
fp1 = open(fileout1, 'w')
fp2 = open(fileout2, 'w')
shorts, longs = read_waves(euvCsv)

fp1.write('\n')
fp1.write('F107 driven model (1) of EUV, using FISM to derive\n')
fp1.write('\n')
fp1.write('model1 = base + Ca * (F107A-50) ^ Pa + Cd * (F107-50) ^ Pd\n')
fp1.write('\n')
fp1.write('bin, short, long, base, Ca, Pa, Cd, Pd, error\n')
fp2.write('\n')
fp1.write('#PARAMS\n')

fp2.write('\n')
fp2.write('F107 driven model (2) of EUV, using FISM to derive\n')
fp2.write('\n')
fp2.write('model2 = base + Ca * (F107A-50) ^ Pa + Cd * (F107-F107A)\n')
fp2.write('\n')
fp2.write('bin, short, long, base, Ca, Pa, Cd, error\n')
fp2.write('\n')
fp2.write('#PARAMS\n')

baseF107 = 50.0
fpert = f107 - baseF107
apert = f107a - baseF107

for iWave in range(0, nWaves):

    if (iWave > 50):
        donorm = False
    else:
        donorm = True
    ir = fismData['irradiances'][:,iWave]

    slopeF, slopeA, inter, pF, pA = fit_model(fpert, apert, ir, donorm = donorm)
    model1 = \
      slopeF * (fpert**pF) + \
      slopeA * (apert**pA) + inter
    model1[model1 < 0.0] = 0.0
    error1 = np.sqrt(np.median(((model1 - ir)/ir*100.0)**2)) 
    #error = calc_error(ir, model1, donorm = donorm)

    sLine = '%d,' % iWave
    sLine = sLine + '%e,' % shorts[iWave]
    sLine = sLine + '%e,' % longs[iWave]
    sLine = sLine + '%e,' % inter
    sLine = sLine + '%e,' % slopeA
    sLine = sLine + '%e,' % pA
    sLine = sLine + '%e,' % slopeF
    sLine = sLine + '%e,' % pF
    sLine = sLine + '%e\n' % error1
    fp1.write(sLine)

    slopeF, slopeA, inter, pF, pA = fit_model2(fpert, apert, ir, donorm = donorm)
    model2 = \
      slopeF * (fpert - apert) + \
      slopeA * (apert**pA) + inter
    model2[model2 < 0.0] = 0.0
    error2 = np.sqrt(np.median(((model2 - ir)/ir*100.0)**2))
    #error = calc_error(ir, model2, donorm = donorm)

    sLine = '%d,' % iWave
    sLine = sLine + '%e,' % shorts[iWave]
    sLine = sLine + '%e,' % longs[iWave]
    sLine = sLine + '%e,' % inter
    sLine = sLine + '%e,' % slopeA
    sLine = sLine + '%e,' % pA
    sLine = sLine + '%e,' % slopeF
    sLine = sLine + '%e\n' % error2
    fp2.write(sLine)
    
    xSize = 7.5
    ySize = 10.0
    dpi = 150
    fig = plt.figure(figsize=(xSize, ySize), dpi = dpi)

    yBot = 0.11
    yTop = 0.05
    yBuf = 0.05
    ax = axes.get_axes_one_column(fig,
                                  2,
                                  yBot,
                                  yTop,
                                  yBuf, xLeft = 0.15)

    ax[0].plot(fismData['times'],  ir, linewidth = 3.0, color = 'k', alpha = 0.3, label = 'FISM')
    ax[0].plot(fismData['times'],  model1, linewidth = 2, color = 'r', label = 'Model1')
    ax[0].plot(fismData['times'],  model2, linewidth = 1, color = 'b', label = 'Model2', alpha = 0.5)
    ax[0].axhline(np.percentile(ir, iLowNorm), color = 'k')
    ax[0].axhline(np.percentile(ir, iLowCutoff), color = 'k', linestyle = ':')
    ax[0].axhline(np.percentile(ir, iHighNorm), color = 'k')
    ax[0].axhline(np.percentile(ir, iHighCutoff), color = 'k', linestyle = ':')
    ax[0].axhline(np.median(ir), color = 'k', linestyle = '--')
    ax[0].set_ylabel('Irradiance (W/m2)')
    ax[0].legend()
    ax[0].set_title('Wavelength Bin : %5.1f - %5.1f nm' % (shorts[iWave], longs[iWave]))
    
    ax[1].scatter(model1, ir, color = 'r', s = 20, label = 'Model1 (%5.1f)' % error1)
    ax[1].scatter(model2, ir, color = 'b', s = 5, label = 'Model2 (%5.1f)' % error2, alpha = 0.5)
    xs = [np.min(ir), np.max(ir)]
    ax[1].plot(xs, xs)
    ax[1].set_xlabel('F107 Derived Model Irradiance (W/m2)')
    ax[1].set_ylabel('FISM Irradiance (W/m2)')
    ax[1].set_aspect(1.0)
    ax[1].axhline(np.percentile(ir, iLowNorm), color = 'k')
    ax[1].axhline(np.percentile(ir, iLowCutoff), color = 'k', linestyle = ':')
    ax[1].axhline(np.percentile(ir, iHighNorm), color = 'k')
    ax[1].axhline(np.percentile(ir, iHighCutoff), color = 'k', linestyle = ':')
    ax[1].axhline(np.median(ir), color = 'k', linestyle = '--')
    ax[1].legend()
    
    outFile = 'euv_%02d.png' % iWave
    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile)
    plt.close(fig)

fp1.close()
fp2.close()
