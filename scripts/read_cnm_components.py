import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt

from astropy.io              import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable




## Multiple (N) Gaussians + offset. ##
 #
 # params list  v    VLSR
 # params float zr   estimated constant zero offset of the data points.
 # params list  h    the array of N estimated heights of the Gaussians.
 # params list  v0   the array of N estimated centers of the Gaussians.
 # params list  w    the array of N estimated halfwidths of the Gaussians.
 #
 # return 1-D-array  tf  The calculated points.
 #
 # version 01/2017
 # author Nguyen Van Hiep ##
def gfunc(v, zr, h, v0, w):
	dp600 = np.float64(0.60056120)
	if(np.isscalar(v)):
		v  = np.array([v], dtype='float64')
	if(np.isscalar(h)):
		h  = np.array([h], dtype='float64')
	if(np.isscalar(v0)):
		v0 = np.array([v0], dtype='float64')
	if(np.isscalar(w)):
		w  = np.array([w], dtype='float64')

	#DETERMINE NR OF GAUSSIANS...
	ng = len(h)
	
	tf = 0.*v + zr
	for i in range(ng):
		if (w[i] > 0.):
			tf = tf + h[i]*np.exp(- ( (v-v0[i])/(dp600*w[i]))**2)

	return tf







#---------- MAIN ----------#

data_path = 'HI_spectra_and_fitted_parameters/HI_spectra/'
src       = '4C+07.13'
data_file = data_path + src + '.fits'

params_path = 'HI_spectra_and_fitted_parameters/fitted_parameters/'
CNM_params  = params_path + 'CNM_fitted_parameters/CNM_Gaussian_parameters.fits'
WNM_params  = params_path + 'WNM_fitted_parameters/WNM_Gaussian_parameters.fits'


## READ CNM components for all sources
data, header = fits.getdata(CNM_params, header=True)
for x in header:
    print (x, header[x])

hdr          = fits.getheader(CNM_params)
for x in hdr:
    print (x, hdr[x])

srcname = data['NAME']
idx     = np.where( srcname == src)[0]
ghgtc   = data['TAU0'][idx]
gcenc   = data['VCEN'][idx]
gwidc   = data['FWHM0'][idx]



## READ SPECTRUM OF A SOURCE
data, header = fits.getdata(data_file, header=True)
hdr          = fits.getheader(data_file)

vlsr         = data['vlsr_tau']
emt          = data['emt']
sigemt       = data['sigemt']

xfitcnm      = data['vlsr_fitcnm']
gfitcnm      = data['fitcnm']
gresid1      = emt - gfitcnm




## Plot Figure
fts  = 36
lbsz = 16
lgds = 14

xlim     = [-65., 40.]
ylime    = [-3., 63.]
ylimeres = [-3., 3.5]
ylima    = [0.3, 1.03]
ylimares = [-0.028, 0.025]

plt.rc('font', weight='bold')
# plt.rc('text', usetex=True)    # uncomment this to make better figure, need to install latex for python
plt.rc('xtick', labelsize=15)
# plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath'] # uncomment this to make better figure, need to install latex for python
mpl.rcParams['axes.linewidth']      = 2.

## Figures
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14,10), sharey=False, sharex=False)

## Absorption
major_xticks = np.arange(-100., 100., 25.)
minor_xticks = np.arange(-100., 100., 5.)
major_yticks = np.arange(0.2, 1.2, 0.1)
minor_yticks = np.arange(0.2, 1.2, 0.05)

ax.plot(vlsr, emt, 'k-', label=r'$e^{-\tau}$')
ax.plot(xfitcnm, gfitcnm, 'r-', label='$Total\ fit$')

divider = make_axes_locatable(ax)
axres   = divider.append_axes('top', size='25%', pad=0, sharex=ax)
ax.figure.add_axes(axres)

for ceni,hgti,widi in zip(gcenc, ghgtc, gwidc):
    ax.plot( vlsr, np.exp(-gfunc(vlsr, 0., [hgti], [ceni], [widi]) ), 'k:', label='')

ax.plot(-100., 1., color='k', ls=':', lw=1, label='$CNM\ components$')
# ax.fill_between(vlsr, emt-sigemt, emt+sigemt, color='gray', label='')

ax.set_ylabel(r'$e^{-\tau}$', fontsize = 18)
ax.set_xlabel(r'$VLSR\ [km\ s^{-1}]$', fontsize = 16)

ax.set_yticks(major_yticks)
ax.set_yticks(minor_yticks, minor=True)
ax.set_xticks(major_xticks)
ax.set_xticks(minor_xticks, minor=True)
ax.tick_params(axis='x', labelsize=lbsz, pad=8)
ax.tick_params(axis='y', labelsize=lbsz)
ax.tick_params(which='both', width=2)
ax.tick_params(which='major', length=6)
ax.tick_params(which='minor', length=3)
ax.legend(loc='lower left', fontsize=lgds)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x',  which='both', bottom=True, top=False, labeltop=False)
ax.axhline(y=ylima[1], xmin=-100, xmax=100, c='k', ls='-.', linewidth=2)

ax.set_xlim(xlim)
ax.set_ylim(ylima)

# Axe - residuals
major_yticks = np.arange(-0.06, 0.06, 0.02)
minor_yticks = np.arange(-0.06, 0.06, 0.01)
axres.plot(vlsr, -sigemt, color='k', ls=':', lw=1, label='')
axres.plot(vlsr, sigemt, color='k', ls=':', lw=1, label='')
axres.plot(xfitcnm, gresid1, color='k', ls='-', label='$Residuals$')

axres.set_ylabel(r'$\mathrm{Residuals}$', fontsize = 16)
axres.set_yticks(major_yticks)
axres.set_yticks(minor_yticks, minor=True)
axres.set_xticks(major_xticks)
axres.set_xticks(minor_xticks, minor=True)
axres.tick_params(axis='x',  which='both', bottom=False, top=True, labelbottom=False)
axres.tick_params(which='both', width=2)
axres.tick_params(which='major', length=6)
axres.tick_params(which='minor', length=3)
axres.grid(False)
axres.spines['bottom'].set_visible(False)
axres.set_xlim(xlim)
axres.set_ylim(ylimares)
plt.savefig(src+'.png', bbox_inches='tight', pad_inches=0.09, format='png', dpi=60)
plt.show()