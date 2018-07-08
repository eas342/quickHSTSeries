import glob
from astropy.io import fits
from astropy.time import Time
import numpy as np
from astropy.table import Table
import pdb

fltList = glob.glob('HST/ob*/*flt.fits')

mjdExp = []
totF, errF = [], []

saveParams = ['PR_INV_L','PR_INV_F','PROPOSID','TDATEOBS','PHOTMODE']
dictOfParamArrays = {}
for oneParam in saveParams:
    dictOfParamArrays[oneParam] = []

for fileInd,oneFile in enumerate(fltList):
    HDUList = fits.open(oneFile)
    primParams = {}
    for oneParam in saveParams:
        primParams[oneParam] = HDUList['PRIMARY'].header[oneParam]
    
    for ind,oneExten in enumerate(HDUList):
        if oneExten.name == 'SCI':
            head = oneExten.header
            expMid = np.mean([head['EXPSTART'],head['EXPEND']])
            mjdExp.append(expMid)
            totF.append(np.sum(oneExten.data))
            if ind + 1 < (len(HDUList) - 1):
                err2D = HDUList[ind + 1].data
                errF.append(np.sqrt(np.sum(err2D**2)))
            else:
                errF.append(np.nan)
            for oneParam in saveParams:
                dictOfParamArrays[oneParam].append(primParams[oneParam])
    HDUList.close()

    if np.mod(fileInd,5) == 0:
        print("Done with file {} out of {}".format(fileInd,len(fltList)))


t = Table()
t['MJD'] = mjdExp
t['Flux'] = totF
t['Flux Err'] = errF

for oneParam in saveParams:
    t[oneParam] = dictOfParamArrays[oneParam]

t.write('flux_totals.csv',overwrite=True)
