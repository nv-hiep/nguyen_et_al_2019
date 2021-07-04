import sys, os
sys.path.insert(0, os.getenv("HOME")+'/Phd@MQ/projects/Dark') # add folder of Class

from common.myimport import *



####################################
#============== MAIN ==============#
####################################

## Path to project
proj_path    = const.HJ2PATH

## Read fit params of HI components from SPONG21 ##
datafile     = const.HJ2PATH + 'data/sponge21/SPONGE_DR1_absparams.fits'
data, header = fits.getdata(datafile, header=True)
ts_sp1       = data['ts']
tkmaxc       = data['tkmax']

for xh in header:
	print xh, header[xh]

# tkmaxc       = 21.866*widcnm**2
Ms           = 4.2*(tkmaxc/ts_sp1 - 1.)

idx = np.where( (Ms>0.) & (Ms<10000.) ) [0]
y   = Ms[idx]
y   = np.sqrt(y)

print 'Med:', np.median(y)

plt.hist(y, 35, density=False, color='k', histtype='step', alpha=1., ls='-', lw=2, label='')
plt.show()

## Filename to save
ret          = {}
ret['SP1_Ms'] = y
filename     = const.HIDIR + 'SPONGE/SP1_Mach_number.npy'
np.save(filename, ret)

sys.exit()