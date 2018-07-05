import QUICFit.QUICFit import QCF

############# Prepping the data #######################

spectra_DIR = './spectra_directory'
save_DIR = './save_directory'

lines_list = np.array([1240,1260,1268,1303,1335,1548])
lines_list_err = np.array([[5,7],[5,2],[3,3],[3,5],[10,10],[5,5]])

QSO_redshift = 5.2

# Spectra should included in Lambda (Angstroms), Flux, Err three-column ascii format
# They should be already power-law normalized for better performance!

spectra = np.loadtxt('./spectra_txt/' + file_spectra, skiprows=1)

################## Running the fitter ################



fitter = QCF.QSOContFitter(spectra[:,0], spectra[:,1], spectra[:,2],
				QSO_lines_list = lines_list, QSO_lines_list_err = lines_list_err,
				redshift_QSO=z_QSO, binning = binning)

fitter.fit(bounds = [1240,1700], show_plots = True,  bin_width = 20, n_sigma = 1.5,
           save_fit_plot = True, directory = './', filename  =  file_spectra[0:10]  +'_fit.pdf')

fitter.save_output(directory = './',
                   filename = file_spectra[0:10]  +'_with_continuum.txt')
