QUICFit
======
QUasar Intrinsic Continuum Fit - A quick, light-weight, automatic spline fitting of the QSO continuum redwards of the Lyman-alpha emission

### Installing QUICFit

Provided you have a working build Python 3.6 with distutils, matplotlib, scipy and numpy, installation should be as simple as

```
$ git clone https://github.com/rameyer/QUICFit
$ cd QUICFit
$ python3 setup.py install
$ cd ..
$ python3
>>> import QUICFit
```

If this fails, you may need more detailed installation instruction. Please contact me! Alternatively, the code is relatively simple and can be used as a local script.

### Getting started

QUICFit is based (for now) on an unique module QUICFit.QSOContFitter, which should be called as 

```
$ python3
>>> import SEDmaker.QSOContFitter as QCF
>>> spectra = np.loadtxt('./DIRECTORY/YOUR_NICE_QSO_SPECTRA.txt')
>>> fitter = QCF(wave = spectra[:,0], flux = spectra[:,1], err = spectra[:,2], redshift_QSO = 5)
```

The observed wavelength (in Angstroms), flux and error array (Arbitrary units) should be supplied at initialization. Other parameters can be supplied as: 
```
	wave: A N array representing the wavelength of the flux datapoints
	flux: A N array with the flux values (power-law normalized)
	std: A N with the std/error of the instrument
	redshift_QSO : The known redshift of the QSO
	binning: Re-binning of the spectra (>= 1)
	QSO_lines_list: The list of expected rest-frame QSO broad lines
	QSO_lines_list_err: The tolerance in rest-frame A of the lines location
```

Once the fitter is initialized, one can run the fit and save the result in a text file containing four columns: wavelength (rest-frame), flux, err, continuum:

```
>>> fitter.fit(bounds = [1240,1700], show_plots = True, chi_min = 0.5, chi_max = 5, 
		   N_steps = 50, bin_width = 20, n_sigma = 1.5, save_fit_plot = True,
		   directory = './save_QSOcont/', filename  =  'QSO_fit.pdf')
>>> fitter.save_output(directory = './save_QSOcont/',filename = 'QSO_continuum.txt')
```

The important parameters here are:

```
	bounds: 1x2 array, the bounds in rest-frame wavelenght where the fit is to be performed
	show_plots: Boolean, whether or not you would like to see the intermediate steps of the fit
	chi_min: Float, minimal chi square for the last spline fit
	chi_max: Float, maximal chi square for the last spline fit
	N_steps: Float, Number of iterations to find the best-fit spline solution
	bin_width: Float, Pixel width of windows to compute the local variance of the pixel-to-pixel flux difference
	n_sigma: Float , Number of sigma from the mean of the distribution of the local flux difference variance to retain continuum pixels
	save_fit_plot: Boolean, whether or not you would like to save a pdf of the final fit complete with observed flux, fit, and residuals
	directory: str , path to the directory where the fit is to be saved
	filename: str, name of the fit pdf file
```


### Improving QUICFit

If you find any error/bug or have suggestions/improvements you'd like to see implemented, please raise an issue or contact me!