import numpy             as np
import matplotlib.pyplot as plt

from   astropy.io        import fits


data_path = 'HI_spectra_and_fitted_parameters/HI_spectra/'
src       = '4C+07.13'
file_path = data_path + src + '.fits'

data, header = fits.getdata(file_path, header=True)
for x in header:
    print (x, header[x])

hdr = fits.getheader(file_path)
for x in hdr:
    print (x, hdr[x])

plt.title( src )
plt.xlabel( 'VLSR [km/s]' )
plt.ylabel( 'e^(-tau)' )

plt.plot(data['vlsr_tau'], data['emt'])
plt.show()