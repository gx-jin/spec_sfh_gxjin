import numpy as np 
import bagpipes as pipes

from astropy.io import fits
# import matplotlib.pyplot as plt

wave_all = np.arange(3000, 10500, 0.7)

hdu = fits.open('/mnt/e/OneDrive/LoTSS/spec_stack_control.fits')
# hdu.info()
spec_control_cen = hdu[1].data['SPEC_CEN_CONTROL']
spec_control_cen_err = hdu[1].data['SPEC_CEN_ERR_CONTROL']

id_control = hdu[1].data['ID_JIN']
hdu1 = fits.open('/mnt/e/OneDrive/LoTSS/spec_stack_rest.fits')

spec_agn_cen = hdu1[1].data['SPEC_CEN'] * 1e-17

spec_agn_cen_err = hdu1[1].data['SPEC_CEN_ERR'] * 1e-17
id_agn = hdu1[1].data['ID_JIN'] 

hdu.close()
hdu1.close()

def bin(spectrum, binn):
    """ Bins up two or three column spectral data by a specified factor. """

    binn = int(binn)
    nbins = len(spectrum) // binn
    binspec = np.zeros((nbins, spectrum.shape[1]))

    for i in range(binspec.shape[0]):
        spec_slice = spectrum[i*binn:(i+1)*binn, :]
        binspec[i, 0] = np.mean(spec_slice[:, 0])
        binspec[i, 1] = np.mean(spec_slice[:, 1])

        if spectrum.shape[1] == 3:
            binspec[i,2] = (1./float(binn)
                            *np.sqrt(np.sum(spec_slice[:, 2]**2)))

    return binspec

def load_spec(ID):
    id0 = int(ID)
    spectrum = np.c_[wave_all,
                     spec_agn_cen[id_agn==id_control[id0]][0],
                     spec_agn_cen_err[id_agn==id_control[id0]][0]]
    mask = (wave_all > 3680) & (wave_all < 9100)
    return bin(spectrum[mask], 1)

galaxy = pipes.galaxy('0', load_spec, photometry_exists=False )

# fig = galaxy.plot()

dblplaw = {}                        
dblplaw["tau"] = (0., 15.)            
dblplaw["alpha"] = (0.01, 1000.)
dblplaw["beta"] = (0.01, 1000.)
dblplaw["alpha_prior"] = "log_10"
dblplaw["beta_prior"] = "log_10"
dblplaw["massformed"] = (7., 15.)
dblplaw["metallicity"] = (0.1, 2.)
dblplaw["metallicity_prior"] = "log_10"

nebular = {}
nebular["logU"] = -3.

dust = {}
dust["type"] = "CF00"
dust["eta"] = 2.
dust["Av"] = (0., 2.0)
dust["n"] = (0.3, 2.5)
dust["n_prior"] = "Gaussian"
dust["n_prior_mu"] = 0.7
dust["n_prior_sigma"] = 0.3

fit_instructions = {}
fit_instructions["redshift"] = (-0.01, 0.01)
fit_instructions["t_bc"] = 0.01
fit_instructions["redshift_prior"] = "Gaussian"
fit_instructions["redshift_prior_mu"] = 0.0
fit_instructions["redshift_prior_sigma"] = 0.001
# fit_instructions["dblplaw"] = dblplaw 
fit_instructions["nebular"] = nebular
fit_instructions["dust"] = dust