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
src       = '3C138'
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




## READ WNM components for all sources
data, header = fits.getdata(WNM_params, header=True)
for x in header:
    print (x, header[x])

hdr          = fits.getheader(WNM_params)
for x in hdr:
    print (x, hdr[x])

srcname = data['NAME']
idx     = np.where( srcname == src)[0]
ghgtw   = data['HGT'][idx]
gcenw   = data['VCEN'][idx]
gwidw   = data['FWHM0'][idx]
fwnm    = data['FWNM'][idx]

## READ SPECTRUM OF A SOURCE
data, header = fits.getdata(data_file, header=True)
hdr          = fits.getheader(data_file)

vlsr         = data['vlsr_TB']
te           = data['TB']
sigtexp      = data['sigTB']

gtb_tot_fit = data['TB_tot_fit']
gtb_wnm_tot = data['TB_WNM_fit']
gtb_cnm_tot = data['TB_CNM_fit']
gresid2     = te - gtb_tot_fit




## Plot
fts  = 36
lbsz = 16
lgds = 14

xlim     = [-65., 40.]
ylime    = [-3., 63.]
ylimeres = [-3., 3.5]
ylima    = [0.3, 1.03]
ylimares = [-0.028, 0.025]

plt.rc('font', weight='bold')
# plt.rc('text', usetex=True)                                                 # uncomment this to make better figure, need to install latex for python
plt.rc('xtick', labelsize=15)
# plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']    # uncomment this to make better figure, need to install latex for python
mpl.rcParams['axes.linewidth']      = 2.

## Figures
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14,10), sharey=False, sharex=False)

## Emission
major_xticks = np.arange(-100., 100., 25.)
minor_xticks = np.arange(-100., 100., 5.)
major_yticks = np.arange(-20., 100., 10.)
minor_yticks = np.arange(-20., 100., 5.)

divider = make_axes_locatable(ax)
axres   = divider.append_axes('bottom', size='25%', pad=0, sharex=ax)
ax.figure.add_axes(axres)
# axres.spines['bottom'].set_visible(False)

ax.plot(vlsr, te, 'k-', lw=1.2, label=r'$T_{exp}$')
ax.plot(vlsr, gtb_tot_fit, 'r-', label=r'$Total\ fit$')
ax.plot(vlsr, gtb_wnm_tot, color='violet', ls='-', zorder=-1, label=r'$Total\ WNM\ fit$')
ax.plot(vlsr, gtb_cnm_tot, color='b', ls='-', zorder=-1, label=r'$Total\ CNM\ fit$')
# ax.fill_between(vlsr, te-sigtexp, te+sigtexp, color='gray', label=r'')
ax.plot(-100., 0., color='k', ls=':', lw=1, label='$WNM\ components$')

for ceni, hgti, widi in zip(gcenw, ghgtw, gwidw):
    ax.plot( vlsr, gfunc(vlsr, 0., [hgti], [ceni], [widi]), 'k:', label='')

ax.set_title('$'+src+': Gaussian$', fontsize=16)
ax.set_ylabel(r'$\mathrm{T_{exp}\ [K]}$', fontsize = 16)

ax.set_yticks(major_yticks)
ax.set_yticks(minor_yticks, minor=True)
ax.set_xticks(major_xticks)
ax.set_xticks(minor_xticks, minor=True)
ax.tick_params(which='y', width=2)
ax.tick_params(which='major', length=6)
ax.tick_params(which='minor', length=3)
# ax.tick_params(axis='x', labelsize=lbsz, pad=8, bottom=False, labelbottom=False)
ax.tick_params(axis='y', labelsize=lbsz)
ax.tick_params(axis='x',  which='both', bottom=False, top=True, labelbottom=False)
ax.grid(False)
ax.set_xticks([])
ax.spines['bottom'].set_visible(False)

ax.set_xlim(xlim)
ax.set_ylim(ylime)
ax.legend(loc='upper left', fontsize=lgds)


# Axe - residuals
major_xticks = np.arange(-100., 100., 25.)
minor_xticks = np.arange(-100., 100., 5.)
major_yticks = np.arange(-6., 6., 2.)
minor_yticks = np.arange(-6., 6., 1.)
axres.plot(vlsr, -sigtexp, color='k', ls=':', lw=1, label='')
axres.plot(vlsr, sigtexp, color='k', ls=':', lw=1, label='')
axres.plot(vlsr, gresid2, color='k', ls='-', label='$Residuals$')

axres.set_ylabel(r'$\mathrm{Residuals}$', fontsize = 16)
axres.set_yticks(major_yticks)
axres.set_yticks(minor_yticks, minor=True)
axres.set_xticks(major_xticks)
axres.set_xticks(minor_xticks, minor=True)
axres.tick_params(axis='x',  which='both', bottom=True, top=False, labelbottom=False)
axres.tick_params(which='both', width=2)
axres.tick_params(which='major', length=6)
axres.tick_params(which='minor', length=3)
axres.grid(False)
axres.spines['top'].set_visible(False) #axres.spines['top'].set_linewidth(0.5)
axres.axhline(y=ylimeres[1], xmin=-100, xmax=100, c='k', ls='-.', linewidth=2)
axres.set_xlim(xlim)
axres.set_ylim(ylimeres)
plt.savefig(src+'_Texp.png', bbox_inches='tight', pad_inches=0.09, format='png', dpi=60)
plt.show()