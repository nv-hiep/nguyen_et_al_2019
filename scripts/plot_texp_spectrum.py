## READ WNM components for all sources
data, header = fits.getdata('spectra_fits/wnm/WNM_params.fits', header=True)
for x in header:
    print x, header[x]

hdr          = fits.getheader('spectra_fits/wnm/WNM_params.fits')
for x in hdr:
    print x, hdr[x]

srcname = data['NAME']
idx     = np.where( srcname == '3C138')[0]
ghgtw   = data['HGT'][idx]
gcenw   = data['VCEN'][idx]
gwidw   = data['FWHM0'][idx]
fwnm    = data['FWNM'][idx]

## READ SPECTRUM OF A SOURCE
data, header = fits.getdata('spectra_fits/3C138.fits', header=True)
hdr          = fits.getheader('spectra_fits/3C138.fits')

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
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
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
    ax.plot( vlsr, cfit.gfunc(vlsr, 0., [hgti], [ceni], [widi]), 'k:', label='')

ax.set_title('$3C138: Gaussian$', fontsize=16)
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
plt.savefig('3C138_Texp.png', bbox_inches='tight', pad_inches=0.09, format='png', dpi=60)
plt.show()