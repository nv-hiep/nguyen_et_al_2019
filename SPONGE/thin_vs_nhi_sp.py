import sys, os
sys.path.insert(0, os.getenv("HOME")+'/Phd@MQ/projects/Dark') # add folder of Class

from common.myimport import *



####################################
#============== MAIN ==============#
####################################

## Infor
infor     = txtdat.read_GNOME_infor_allsrc(fname =  const.HJ2PATH+'source/idlCarl/data/80GNOMES_src_basic_infor.csv', asarray=True)
okyn      = infor['okyn']
srclist   = infor['src_cld']

# idx       = np.where( okyn==1 ) [0]
idx       = np.where( (okyn==1) & (srclist != 'SRC14_T') & (srclist != 'SRC29_R') ) [0]
srclist   = srclist[idx]
allsrc    = infor['src'][idx]
clds      = infor['cloud'][idx]
xv1       = infor['v1'][idx]
xv2       = infor['v2'][idx]
xv3       = infor['v3'][idx]
xv4       = infor['v4'][idx]
xl        = infor['glong'][idx]
xb        = infor['glat'][idx]
ra        = infor['ra'][idx]
dec       = infor['dec'][idx]
okyn      = infor['okyn'][idx]





## Path to project
# cols = ['src','l', 'b', 'cnm','cnm_er','wnm','wnm_er','nhi','nhi_er','thin','thin_er']
# fmt  = ['s',  'f', 'f', 'f',    'f',   'f',  'f',     'f',   'f',     'f',    'f'    ]
dat     = txtdat.read_info_sponge_30src(fname = const.HJ2PATH + 'data/sponge21/30src_claire.txt', asarray=True)

spsrc   = dat['src'] 
xlsp    = dat['l']
xbsp    = dat['b']
nhi     = dat['nhi']
signhi  = dat['nhi_er']
thin    = dat['thin']
sigthin = dat['thin_er']

ratio   = nhi/thin
sigrat  = stats.uncertainty_of_ratio(nhi, thin, signhi, sigthin)

plt.plot(thin, ratio, 'k.')
for i,sc in enumerate(spsrc):
	plt.annotate('('+sc+')', xy=(thin[i], ratio[i]), xycoords='data', xytext=(-50.,30.), textcoords='offset points', arrowprops=dict(arrowstyle="->"),fontsize=14)

plt.show()


## PLOT - ratio vs thin  for all 30  SPONGE soures ##
fts  = 16
lbsz = 14
lgds = 10

plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
mpl.rcParams['axes.linewidth']      = 2.


fig, ax = plt.subplots( 1, 1, figsize=(10,10) )
# axes     = axs.flatten()

# ax = axs[0]


major_xticks = np.arange(160., 215., 10.)
minor_xticks = np.arange(160., 215., 5.)
major_yticks = np.arange(-40., 30., 10.)
minor_yticks = np.arange(-40., 30., 5.)

ax.plot(xlsp, xbsp, 'kx')
ax.plot(xl, xb, 'r+')
ax.plot([160., 215.], [0., 0.], zorder=-1, color='grey', ls=':', lw=2)

## Plot regions
ax.plot([216.8, 200.7, 187.8, 159.2, 159.2, 180., 202.39, 206.2, 216.8], [7.9, 20.9, 2.8, 2.8, -38.38, -38.38, -10.36, -12.4, 7.9], color='grey', ls='--', zorder=-1)

## Plot clouds
print 'plot: Taurus'
glx = [162., 180.]
gbx = [-22., -10.]
ax.plot( [glx[0], glx[1]], [gbx[1], gbx[1]], 'k--', zorder=-1, label='')
ax.plot( [glx[0], glx[0]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[1], glx[1]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[0], glx[1]], [gbx[0], gbx[0]], 'k--', zorder=-1)
ax.text( 173.9, -21.52, '$Taurus$' )

print 'plot: NGC 2264'
glx = [200., 204.]
gbx = [0., 4.]
ax.plot( [glx[0], glx[1]], [gbx[1], gbx[1]], 'k--', zorder=-1, label='')
ax.plot( [glx[0], glx[0]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[1], glx[1]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[0], glx[1]], [gbx[0], gbx[0]], 'k--', zorder=-1)
ax.text( 203.3, 4.08, '$NGC2264$' )

print 'plot: Rosette'
glx = [205., 209.]
gbx = [-3., -1.]
ax.plot( [glx[0], glx[1]], [gbx[1], gbx[1]], 'k--', zorder=-1, label='')
ax.plot( [glx[0], glx[0]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[1], glx[1]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[0], glx[1]], [gbx[0], gbx[0]], 'k--', zorder=-1)
ax.text( 212.7, -2., '$Rosette$' )

print 'plot: California'
glx = [154., 164.]
gbx = [-15., -5.]
ax.plot( [glx[0], glx[1]], [gbx[1], gbx[1]], 'k--', zorder=-1, label='')
ax.plot( [glx[0], glx[0]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[1], glx[1]], [gbx[0], gbx[1]], 'k--', zorder=-1)
ax.plot( [glx[0], glx[1]], [gbx[0], gbx[0]], 'k--', zorder=-1)
ax.text( 165.6, -4.4, '$California$' )
## ENd - Clouds

ax.set_ylabel(r'$\mathrm{Galactic\ latitude\ [^{o}]}$', fontsize=fts, fontweight='normal')
ax.set_xlabel(r'$\mathrm{Galactic\ longitude\ [^{o}]}$', fontsize=fts, fontweight='normal')

ax.set_yticks(major_yticks)
ax.set_yticks(minor_yticks, minor=True)

ax.set_xticks(major_xticks)
ax.set_xticks(minor_xticks, minor=True)

ax.tick_params(axis='x', labelsize=lbsz, pad=4)
ax.tick_params(axis='y', labelsize=lbsz)
ax.tick_params(which='both', width=2)
ax.tick_params(which='major', length=5)
ax.tick_params(which='minor', length=2)
ax.grid(False)

ax.set_xlim(215., 159.)
ax.set_ylim(-42., 20.)

for i,sc in enumerate(spsrc):
	ax.annotate('(' + sc + ')', xy=(xlsp[i], xbsp[i]), xycoords='data', xytext=(-50.,30.), textcoords='offset points', arrowprops=dict(arrowstyle="->"),fontsize=14)

## Legend
# plt.legend(loc='lower left', fontsize=lgds, handletextpad=-0.2, numpoints=1)
# plt.legend(loc='lower left', fontsize=lgds)


# plt.savefig(const.PLOTDIR+'fcnm_vs_nhi.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=80)
plt.show()

slist   = ['3C111', '3C123', '3C133']
srclist = spsrc.tolist()
idx     = [ srclist.index(x) for x in slist]


ret            = {}
ret['src']     = spsrc[idx]
ret['thin']    = thin[idx]
ret['sigthin'] = sigthin[idx]
ret['nhi']     = nhi[idx]
ret['signhi']  = signhi[idx]

ret['ratio']   = ratio[idx]
ret['sigrat']  = sigrat[idx]

ret['l']       = xlsp[idx]
ret['b']       = xbsp[idx]
ret['survey']  = [1]*len(slist)

filename       = const.HIDIR+'read_results/both_GV/model_selection/ratio_vs_thin_sponge.npy'
np.save(filename, ret)


idx = np.where( ret['b'] < -5. )[0]
print idx
print ret['src'][idx]


ret['src']     = ret['src'][idx]
ret['thin']    = ret['thin'][idx]
ret['sigthin'] = ret['sigthin'][idx]
ret['nhi']     = ret['nhi'][idx]
ret['signhi']  = ret['signhi'][idx]

ret['ratio']   = ret['ratio'][idx]
ret['sigrat']  = ret['sigrat'][idx]

ret['l']       = ret['l'][idx]
ret['b']       = ret['b'][idx]
ret['survey']  = [1]*len(idx)

print len(ret['src'])

filename       = const.HIDIR+'read_results/both_GV/model_selection/ratio_vs_thin_sponge_Taurus.npy'  ## 3 srcs ['3C111', '3C123', '3C133']
np.save(filename, ret)