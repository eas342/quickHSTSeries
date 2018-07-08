import glob
from astropy.io import fits
from astropy.time import Time
import numpy as np
from astropy.table import Table
import pdb

fltList = glob.glob('HST/ob*/*flt.fits')

mjdExp = []
totF, errF = [], []

for fileInd,oneFile in enumerate(fltList):
    HDUList = fits.open(oneFile)
    for ind,oneExten in enumerate(HDUList):
        if oneExten.name == 'SCI':
            head = oneExten.header
            expMid = np.mean([head['EXPSTART'],head['EXPEND']])
            mjdExp.append(expMid)
            pdb.set_trace()
            totF.append(np.sum(oneExten.data))
            if ind + 1 < (len(HDUList) - 1):
                err2D = HDUList[ind + 1].data
                errF.append(np.sqrt(np.sum(err2D**2)))
            else:
                errF.append(np.nan)
    HDUList.close()

    if np.mod(fileInd,5) == 0:
        print("Done with file {} out of {}".format(fileInd,len(fltList)))


t = Table()
t['MJD'] = mjdExp
t['Flux'] = totF
t['Flux Err'] = errF
t.write('flux_totals.csv',overwrite=True)
