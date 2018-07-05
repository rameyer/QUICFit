import numpy as np
import pkg_resources
import QUICFit.QUICFit as QCF

############# Prepping the data #######################

spectra_DIR = './spectra_directory'
save_DIR = './save_directory'

lines_list = np.array([1240,1263,1270,1303,1330,1335,1400,1440,1548])
lines_list_err = np.array([[5,7],[5,2],[3,3],[3,5],[5,2],[2,10],[10,10],[5,5],[5,5]])

QSO_redshift = 5.79

# Spectra should included in Lambda (Angstroms), Flux, Err three-column ascii format
# They should be already power-law normalized for better performance!
# The next line uses an example QSO, but feel free to change any of your choice!
spectra = np.loadtxt(pkg_resources.resource_filename('QUICFit','/data/J0002+2550_ESI_arch.txt'), skiprows=1)

################## Running the fitter ################

fitter = QCF.QSOContFitter(wave=spectra[:,0], flux=spectra[:,1], err=spectra[:,2],
				QSO_lines_list = lines_list, QSO_lines_list_err = lines_list_err,
				redshift_QSO=QSO_redshift, binning = 1)

fitter.fit(bounds = [1240,1700], show_plots = True,  bin_width = 20, n_sigma = 2,
           save_fit_plot = True, directory = './', filename  = 'J0002+2550_fit.pdf')

fitter.save_output(directory = './',
                   filename =  'J0002+2550_with_continuum.txt')
