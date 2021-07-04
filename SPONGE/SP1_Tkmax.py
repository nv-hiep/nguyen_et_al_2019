import sys, os
sys.path.insert(0, os.getenv("HOME")+'/Phd@MQ/projects/Dark') # add folder of Class

from common.myimport import *



####################################
#============== MAIN ==============#
####################################

## Path to project
proj_path    = const.HJ2PATH

## Read fit params of HI components from SPONG21 ##
datafile     = const.HJ2PATH + 'data/sponge21/SPONGE_DR1_emparams.fits'
data, header = fits.getdata(datafile, header=True)

print header

for x in header:
	print x, header[x]

srcnames     = data['names']
# for x in srcnames:
# 	print x
sys.exit()

widwnm       = data['FWHM0']
# tkmaxc       = data['tkmax']

kk = 0
kall = 0
for wwi in widwnm:
	kall += 1
	if(wwi < 3.2):
		print wwi
		kk = kk + 1

print kk
print kall

sys.exit()

tkmaxw       = 21.866*widwnm**2
# Ms           = 4.2*(widwnm/ts_sp1 - 1.)

# idx = np.where( (Ms>0.) & (Ms<10000.) ) [0]
# y   = Ms[idx]
# y   = np.sqrt(y)

# print 'Med:', np.median(y)

plt.hist(tkmaxw, 100, density=False, color='k', histtype='step', alpha=1., ls='-', lw=2, label='')
plt.show()

## Filename to save
ret               = {}
ret['SP1_Tkmaxw'] = tkmaxw
filename          = const.HIDIR + 'SPONGE/SP1_Tkmaxw.npy'
np.save(filename, ret)

sys.exit()