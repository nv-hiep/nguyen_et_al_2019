## READ CNM components for all sources
data, header = fits.getdata('spectra_fits/cnm/CNM_params.fits', header=True)
for x in header:
    print x, header[x]

hdr          = fits.getheader('spectra_fits/cnm/CNM_params.fits')
for x in hdr:
    print x, hdr[x]

srcname = data['NAME']
idx     = np.where( srcname == '3C138')[0]
ghgtc   = data['TAU0'][idx]
gcenc   = data['VCEN'][idx]
gwidc   = data['FWHM0'][idx]



## READ SPECTRUM OF A SOURCE
data, header = fits.getdata('spectra_fits/3C138.fits', header=True)
hdr          = fits.getheader('spectra_fits/3C138.fits')

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
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
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
    ax.plot( vlsr, np.exp(-cfit.gfunc(vlsr, 0., [hgti], [ceni], [widi]) ), 'k:', label='')

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
plt.savefig('3C138.png', bbox_inches='tight', pad_inches=0.09, format='png', dpi=60)
plt.show()