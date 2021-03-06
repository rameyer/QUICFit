import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join


def gaussian(x, Amp, mu, sigma):
	''' A simple gaussian with the usual parameters
	Inputs:
		x : the input values
		Amp: the Amplitude of the gaussian
		mu : The mean
		sigma: The square root of the variance
	Returns:
		g(x) with the above parameters
	'''

	return Amp * np.exp(-(x - mu) ** 2 / (2 * sigma * sigma))


def matched_filter(u, f, e):
	''' Matched Filtering: e.g. Hewett 1985, Bolton 2004. 
	Inputs:
		u: Kernel/filter
		f: flux of the spectra
		e: error array (std of the flux)
	Returns:
		SNR: Signal to Noise array for all flux values 
		(except start and end which are null due to zero padding)
	'''
	width = int(len(u) / 2)

	Cj1 = np.array([np.sum(u * f[j - width:j + width] / e[j - width:j + width] ** 2)
	                for j in range(width, int(len(f) - width))])
	Cj2 = np.array([np.sum(u ** 2 / e[j - width:j + width] ** 2)
	                for j in range(width, int(len(f) - width))])

	SNR = np.zeros(len(f))
	SNR[width: len(f) - width] = Cj1 / np.sqrt(Cj2)

	return SNR


def running_median(array, width):
	''' Perform a running median with a window of size 2*width
		The first and last N(width) values are 0 and should be discarded
	Inputs:
		array: The array on which to run the median
		widht: Half the window width
	Returns:
		median: The runned median
	'''
	median = np.zeros(len(array))

	median[width:len(array) - width] = np.array([np.median(array[j - width:j + width])
	                                             for j in
	                                             range(width, len(array) - width)])

	return median


def running_average(array, width):
	''' Perform a running average with a window of size 2*width
		The first and last N(width) values are 0 and should be discarded
	Inputs:
		array: The array on which to run the average
		widht: Half the window width
	Returns:
		median: The runned average
	'''
	average = np.zeros(len(array))

	average[width:len(array) - width] = np.array([np.mean(array[j - width:j + width])
	                                              for j in
	                                              range(width, len(array) - width)])

	return average


def running_std(array, width):
	''' Perform a running standard deviation with a window of size 2*width
		The first and last N(width) values are 0 and should be discarded
	Inputs:
		array: The array on which to run the std
		widht: Half the window width
	Returns:
		median: The runned std
	'''
	std = np.zeros(len(array))

	std[width:len(array) - width] = np.array([np.std(array[j - width:j + width])
	                                          for j in range(width, len(array) - width)])

	return std


class QSOContFitter(object):
	'''
	The QSO Continuum fitter. Provides an automated pseudo-continuum spline 
	fitting for all continuum and broad emission redwards of the Lyman alpha
	emission. The useinterface provides three key steps: initialization, 
	fitting and eventually saving the output values of the spline.
	The fitting procedure relies on two key aspects: continuum filtering via 
	matched-filter and statistical selection of the continuum, and 
	semi-automated spline fitting where half of the knot points are forced to
	be close the the expected QSO broad lines supplied by the user. 
	
	Minimal example:
	lines_list = np.array([1240,1260,1265,1303,1335,1400,1549])
	lines_err = np.array([[5,5],[5,1],[3,5],[5,5],[5,5],[5,5],[5,5]])
	spectra = np.loadtxt('./QSO_spectra_txt/',skiprows =1)
	
	Fitter = QSOContFitter(spectra[:,0], spectra[:,1], spectra[:,2], 
				QSO_lines_list = lines_list, QSO_lines_list_err = lines_err,
				redshift_QSO=6, binning = 1)
	Fitter.fit(bounds = [1220,1600], show_plots = True, chi_min = 0.5, 
			   chi_max = 5, N_steps = 50, bin_width = 20, n_sigma = 2, 
			   save_fit_plot = True, directory = './save/', 
			   filename  =  file_spectra[0:10]  +'_fit.pdf')

	Fitter.save_output(directory = './save/',filename = 'QSO_continuum.txt')
	'''

	def __init__(self, wave, flux, err, redshift_QSO, binning=1,
	             QSO_lines_list=np.array([1240, 1260, 1303, 1335, 1400, 1549]),
	             QSO_lines_list_err=np.array([[5,5], [5,5], [5,5], [5.5], [5,5], [5,5]])):
		'''
		Initializes the fitter with the input flux and key parameters. 
		Input:
			wave: A N array representing the wavelength of the flux datapoints
			flux: A N array with the flux values (power-law normalized)
			std: A N with the std/error of the instrument
			redshift_QSO : The known redshift of the QSO
			binning: Re-binning of the spectra (>= 1)
			QSO_lines_list: The list of expected rest-frame QSO broad lines
			QSO_lines_list_err: The tolerance in rest-frame A of the lines
								location
		Returns:
		    - : Re-bin the data and stores the relevant parameters.
		'''

		assert len(wave) == len(
			flux), 'Wavelength and flux array have a discrepant length.'
		assert len(err) == len(
			flux), 'Wavelength and error array have a discrepant length.'

		self.original_wave = wave
		self.original_flux = flux
		self.original_std = err

		self.wave = wave
		self.flux = flux
		self.std = err

		assert redshift_QSO > 0, 'Please input a physical value for the QSO redshift.'

		self.redshift_QSO = redshift_QSO

		assert len(QSO_lines_list) == len(QSO_lines_list_err), \
			'QSO broad lines and errors list have discrepant length!'

		self.QSO_lines_list = QSO_lines_list
		self.QSO_lines_err = QSO_lines_list_err
		self.knots = []
		self.err_knots_lower = []
		self.err_knots_higher = []

		assert binning >= 1, 'Please enter a valid binning value (>=1)'
		assert round(binning) == binning, 'Please enter a valid binning value (integer)'
		self.binning = binning
		if self.binning > 1:
			print('Re-binning the data...')
			self._bin_data()
		else:
			print('No re-binning required.')

	def fit(self, bounds=[0, np.inf], show_plots=False, n_sigma=1, bin_width=10,
	        save_fit_plot=False, directory=None, filename=None):

		self._minmax_wavelength(bounds)

		self.continuum_indices = self._find_continuum(gaussian_amp=-1,
		                                              gaussian_width=10, n_sigma=n_sigma,
		                                              bin_width=bin_width,
		                                              show_plots=show_plots)

		self.spl = self._fit_spline_from_knots()

		if show_plots or save_fit_plot:
			self._plot_spline_and_residuals(show_plot=show_plots,
			                                save_plot=save_fit_plot,
			                                directory=directory,
			                                filename=filename)
		return

	def save_output(self, directory, filename, save_binned_data=True):
		if directory is None or filename is None:
			if directory is None:
				print(
					'Directory is not defined, output pseudo-continuum will not be saved.')
			elif filename is None:
				print(
					'Filename not specified, output pseudo-continuum will not be saved.')
			return
		else:
			if save_binned_data:
				np.savetxt(join(directory, filename),
				           np.transpose([self.wave[self._min_index:self._max_index],
				                         self.flux[self._min_index:self._max_index],
				                         self.std[self._min_index:self._max_index],
				                         self.spl(self.wave[
				                                  self._min_index:self._max_index])]))
			else:
				np.savetxt(join(directory, filename),
				           np.transpose([self.original_wave[
				                         self._min_index * self.binning:self._max_index * self.binning],
				                         self.original_flux[
				                         self._min_index * self.binning:self._max_index * self.binning],
				                         self.original_std[
				                         self._min_index * self.binning:self._max_index * self.binning],
				                         self.spl(self.original_wave[
				                                  self._min_index * self.binning:self._max_index * self.binning])]))

	def _minmax_wavelength(self, bounds):

		assert bounds[0] < bounds[
			1], 'Lower bound must be strictly inferior to superior bound.'

		if self.wave[-1] > (1 + self.redshift_QSO) * bounds[0]:
			self._min_index = np.min(
				np.where(self.wave > bounds[0] * (1 + self.redshift_QSO)))
		else:
			self._min_index = 1
			print(
				'Lower bound inferior to lowest wavelength. Taking wavelength array minimum instead.')
		if self._min_index == 0:
			self._min_index = 1

		if self.wave[-1] > (1 + self.redshift_QSO) * bounds[1]:
			self._max_index = np.min(
				np.where(self.wave > bounds[1] * (1 + self.redshift_QSO)))
		else:
			self._max_index = len(self.wave) - 1
			print(
				'Higher bound superior to maximum wavelength. Taking wavelength array maximum instead.')

	def _find_continuum(self, gaussian_amp=-0.5, gaussian_width=20,
	                    n_sigma=1, bin_width=10, show_plots=False):
		# What is a continuum? Average noise-jumps from pixel to pixel

		kernel = gaussian(np.linspace(0, 20, 20), gaussian_amp, 10, gaussian_width)

		SNR = matched_filter(kernel, self.flux[self._min_index:self._max_index]
		                     - running_median(self.flux, 50)[
		                       self._min_index:self._max_index],
		                     self.std[self._min_index:self._max_index])

		snr_indices = np.where(SNR < 3)

		if show_plots == True:
			#plt.plot(self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
			 #        SNR, '--r')
			plt.plot(self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
			         self.flux[self._min_index:self._max_index], '--c',label='All pixels')
			plt.plot(
				self.wave[self._min_index:self._max_index][snr_indices] / (1 + self.redshift_QSO),
				self.flux[self._min_index:self._max_index][snr_indices], 'k',label=r'$S/N_{matched-filter}<3$')
			plt.title('Matched-filtering of large absorptions',fontsize = 16)
			plt.legend(frameon=False,fontsize=14)
			plt.xlabel(r'Rest-frame wavelength [$\AA$]',fontsize = 16)
			plt.ylabel('Arbitrary',fontsize = 16)
			plt.show()

		print(self._min_index, self._max_index)
		self.dflux_per_dpix = self.flux[self._min_index:self._max_index] - \
		                      self.flux[self._min_index - 1:self._max_index - 1]

		self.estimated_std = running_std(self.flux[1::] -
		                                 self.flux[0:-1], bin_width)[
		                     self._min_index:self._max_index]

		self.average_std = running_average(self.std, bin_width)[
		                   self._min_index:self._max_index]

		temp_std_ratio = np.array([est / average for est, average, s, f, e in
		                           zip(self.estimated_std, self.average_std, SNR,
		                               self.flux[self._min_index:self._max_index],
		                               self.std[self._min_index:self._max_index])
		                           if
		                           average > 0 and est > 0 and s < 3 and f / e >= 5 and est / average < 20])

		hist, bin_edges = np.histogram(temp_std_ratio, 200, normed=True)
		bins_centers = bin_edges[0:-1] + (bin_edges[0:-1] - bin_edges[1::]) / 2.

		hist = hist / np.max(hist)
		index_max = np.where(hist == 1)[0][0]

		parameters_gaussian = curve_fit(gaussian, bins_centers[0:int(index_max * 1.2)],
		                                hist[0:int(index_max * 1.2)],
		                                p0=[1, np.max(bins_centers[int(index_max * 1.2)]),
		                                    0.4],
		                                bounds=([0.5, 0, 0], [1.0, np.inf, np.inf]))[0]

		if show_plots:
			plt.figure(figsize=(6, 6))
			plt.plot(bins_centers, gaussian(bins_centers, parameters_gaussian[0],
			                                parameters_gaussian[1],
			                                parameters_gaussian[2], ), '--c',
			         label='Semi-gaussian fit')
			plt.vlines(x=parameters_gaussian[1] + n_sigma * parameters_gaussian[2],
			           ymin=0, ymax=0.8, color='r', linestyle='dashed',
			           label='Continuum cut')
			plt.step(bins_centers, hist / np.max(hist), where='mid', color='k')
			plt.xlim(0, 5)
			plt.ylabel(r'$N_{pix}$ [normalized]', fontsize=20)
			plt.xlabel(r'$\Sigma_\Delta = \sigma_\Delta / \sigma$', fontsize=20)
			plt.legend(fontsize=18, frameon=False)
			plt.xticks(fontsize=18)
			plt.yticks(fontsize=18)
			ax = plt.gca()
			ax.tick_params(direction='in', which='both', labelsize=18)
			plt.show()

		continuum_indices = np.array([temp_k for temp_k in
		                              range(len(
			                              self.wave[self._min_index:self._max_index - 1]))
		                              if (self.estimated_std[temp_k] < (
						parameters_gaussian[1] +
						n_sigma * parameters_gaussian[2]) * self.average_std[temp_k])])

		continuum_indices = np.array([k for k in continuum_indices if SNR[k] < 3 and
		                              self.std[self._min_index:self._max_index][
			                              k] <= 0.3])

		print('Continuum indices found. Proceding to spline fitting.')

		if show_plots == True:
			#plt.plot(self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
			#        self.estimated_std, 'g')
			#plt.plot(self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
			#         (parameters_gaussian[1] + n_sigma * parameters_gaussian[
			#	         2]) * self.average_std, 'k')
			plt.plot(self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
			         self.flux[self._min_index:self._max_index], '--c',
			         label = 'Original spectra')
			plt.plot(self.wave[self._min_index:self._max_index][continuum_indices] / (
						1 +self.redshift_QSO),
			         self.flux[self._min_index:self._max_index][continuum_indices], 'k',
			         label = 'Continuum pixels')
			plt.title('Final continuum pixel selection',fontsize = 16)
			plt.xlabel(r'Rest-frame wavelength [$\AA$]',fontsize = 16)
			plt.ylabel('Flux',fontsize = 16)
			plt.legend(frameon=False, fontsize=14)
			plt.show()

		return continuum_indices

	def _fit_spline_from_knots(self):
		if self.knots == []:
			self._initialize_knots()

		# define a small epsilon to avoid knots to be sitting in neighbouring pixels
		epsilon = 1

		bounds_knots = [[k - low + epsilon, k + high - epsilon] for k, low, high
		                in zip(self.knots, self.err_knots_lower, self.err_knots_higher)]

		res = minimize(self._residuals_knots, x0=self.knots, method='L-BFGS-B',
		               bounds=bounds_knots)

		self.knots = res.x

		print('Fitted rest-frame knots are: ', self.knots / (1 + self.redshift_QSO))

		return self._spline_from_knots(self.knots)

	def _residuals_knots(self, knots):

		return self._spline_from_knots(knots).get_residual()

	def _spline_from_knots(self, knots):

		spline = LSQUnivariateSpline(
			x=self.wave[self._min_index:self._max_index][self.continuum_indices],
			y=self.flux[self._min_index:self._max_index][self.continuum_indices],
			t=knots,
			w=1. / self.std[self._min_index:self._max_index][self.continuum_indices],
			k=3)
		return spline

	def _initialize_knots(self):

		lines_err_cropped = np.array([[line * (1 + self.redshift_QSO),
		                               low_err * (1 + self.redshift_QSO),
		                               high_err * (1 + self.redshift_QSO)]
		                              for line, low_err, high_err in
		                              zip(self.QSO_lines_list, self.QSO_lines_err[:, 0],
		                                  self.QSO_lines_err[:, 1])
		                              if (line * (1 + self.redshift_QSO) < self.wave[
				self._max_index] and
		                                  line * (1 + self.redshift_QSO) > self.wave[
			                                  self._min_index])])

		self.knots = []
		self.err_knots_lower = []
		self.err_knots_higher = []

		for i in range(len(lines_err_cropped)):
			if i == 0:
				# Mid-point between inferior wavelength range bound & 1st QSO line
				self._append_closest_to_knots_err(
					0.5 * (self.wave[self._min_index] + lines_err_cropped[i, 0]),
					0.5 * (lines_err_cropped[i, 0] - self.wave[self._min_index]),
					0.5 * (lines_err_cropped[i, 0] - self.wave[self._min_index]) -
					lines_err_cropped[i, 1])
				# Actual QSO line knot point
				self.knots.append(lines_err_cropped[i, 0])
				self.err_knots_lower.append(lines_err_cropped[i, 1])
				self.err_knots_higher.append(lines_err_cropped[i, 2])
				# Next knot mid-point
				self._append_closest_to_knots_err(
					0.5 * (lines_err_cropped[i + 1, 0] + lines_err_cropped[i, 0]),
					0.5 * (lines_err_cropped[i + 1, 0] - lines_err_cropped[i, 0]) -
					lines_err_cropped[i, 2],
					0.5 * (lines_err_cropped[i + 1, 0] - lines_err_cropped[i, 0]) -
					lines_err_cropped[i + 1, 1])
			elif i == len(lines_err_cropped) - 1:
				# Last QSO line knot
				self.knots.append(lines_err_cropped[i, 0])
				self.err_knots_lower.append(lines_err_cropped[i, 1])
				self.err_knots_higher.append(lines_err_cropped[i, 2])
				# Mid-point between superior wavelength range bound & last QSO line
				self._append_closest_to_knots_err(
					0.5 * (self.wave[self._max_index] + lines_err_cropped[i, 0]),
					0.5 * (self.wave[self._max_index] - lines_err_cropped[i, 0]) -
					lines_err_cropped[i, 2],
					0.5 * (self.wave[self._max_index] - lines_err_cropped[i, 0]))
			else:
				# Actual QSO line knot point
				self.knots.append(lines_err_cropped[i, 0])
				self.err_knots_lower.append(lines_err_cropped[i, 1])
				self.err_knots_higher.append(lines_err_cropped[i, 2])
				# Next knot mid-point
				self._append_closest_to_knots_err(
					0.5 * (lines_err_cropped[i + 1, 0] + lines_err_cropped[i, 0]),
					0.5 * (lines_err_cropped[i + 1, 0] - lines_err_cropped[i, 0]) -
					lines_err_cropped[i, 2],
					0.5 * (lines_err_cropped[i + 1, 0] - lines_err_cropped[i, 0]) -
					lines_err_cropped[i + 1, 1])

		self.knots = np.array(self.knots)
		self.err_knots_higher = np.array(self.err_knots_higher)
		self.err_knots_lower = np.array(self.err_knots_lower)
		print('Initialized knots to ', self.knots)
		print('Wave bounds', self.wave[self._min_index], self.wave[self._max_index])
		print('Err knots lower to   ', self.err_knots_lower)
		print('err knots higher to  ', self.err_knots_higher)

	def _append_closest_to_knots_err(self, wavelength, tol_low, tol_high):
		possible_values = self._find_continuum_values_in_range(wavelength - tol_low,
		                                                       wavelength + tol_high)

		if possible_values != [] and len(possible_values) >= 3:
			n = len(possible_values)
			self.knots.append(possible_values[int(1 + 0.5 * n)])
			self.err_knots_lower.append(
				possible_values[int(1 + 0.5 * n)] - np.min(possible_values))
			self.err_knots_higher.append(
				np.max(possible_values) - possible_values[int(1 + 0.5 * n)])
		else:
			print(
			'Knot with rest-frame wavelength ', wavelength / (1 + self.redshift_QSO),
			' not initialized to respect Schoenberg-Whitney conditions. Check for sufficient continuum coverage in the region.')

	def _find_continuum_values_in_range(self, wavemin, wavemax):
		temp = np.array(
			[l for l in self.wave[self._min_index:self._max_index][self.continuum_indices]
			 if (l < wavemax and l > wavemin)])
		return temp

	def _1Dspline(self, chi2_spline):

		spline = UnivariateSpline(
			self.wave[self._min_index:self._max_index][self.continuum_indices],
			self.flux[self._min_index:self._max_index][self.continuum_indices],
			w=1. / self.std[self._min_index:self._max_index][self.continuum_indices],
			k=3, s=len(self.continuum_indices) * chi2_spline)

		d_spline = spline.derivative(n=1)(self.wave)[self._min_index:self._max_index]

		maximas_indices = np.array(
			[k for k in range(len(self.wave[self._min_index:self._max_index]) - 1)
			 if (d_spline[k] > 0 and d_spline[k + 1] < 0)])

		maximas_indices = np.int_(np.concatenate([maximas_indices,
		                                          argrelmax((d_spline < 0) * d_spline)[
			                                          0]]))
		return spline, self.wave[self._min_index:self._max_index][maximas_indices]

	def _plot_spline_and_residuals(self, show_plot=True, save_plot=False, directory=None,
	                               filename=None):

		fig = plt.figure(figsize=(16, 5))
		ax = plt.subplot(2, 1, 1)

		ax.plot(self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
		        self.flux[self._min_index:self._max_index], 'k',
		        label='Power-law normalized QSO spectra')
		ax.plot(self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
		        self.spl(self.wave[self._min_index:self._max_index]), '--c',
		        label='Automated Spline Fit', linewidth=2)
		ax.set_ylim(np.max([0.0, np.min(self.flux[self._min_index:self._max_index])]),
		            np.min([3.0, np.max(self.flux[self._min_index:self._max_index])]))
		ax.set_xlim(self.wave[self._min_index] / (1 + self.redshift_QSO),
		            self.wave[self._max_index] / (1 + self.redshift_QSO))
		ax.set_ylabel('Normalized Flux', fontsize=18)
		ax.tick_params(direction='in')
		ax.xaxis.set_ticklabels([])
		plt.ylim(0.5, 1.8)
		ax.yaxis.set_label_coords(-0.06, 0.5)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		ax.legend(fontsize=14, loc=2, frameon=False, ncol=2)
		ax = plt.gca()
		ax.tick_params(direction='in', which='both', labelsize=16)
		ax = plt.subplot(2, 1, 2)
		residuals = (self.flux[self._min_index:self._max_index][self.continuum_indices]
		             - self.spl(
					self.wave[self._min_index:self._max_index][self.continuum_indices]))
		ax.scatter(self.wave[self._min_index:self._max_index][self.continuum_indices] / (
					1 + self.redshift_QSO),
		           residuals, color='k', s=5)
		ax.set_ylim(-np.max(np.abs(residuals)), np.max(np.abs(residuals)))
		ax.set_xlim(self.wave[self._min_index] / (1 + self.redshift_QSO),
		            self.wave[self._max_index] / (1 + self.redshift_QSO))
		ax.fill_between(
			self.wave[self._min_index:self._max_index] / (1 + self.redshift_QSO),
			2 * self.std[self._min_index:self._max_index],
			-2 * self.std[self._min_index:self._max_index],
			color='r', alpha=0.5)
		ax.tick_params(direction='in')
		ax.yaxis.set_ticks([-0.25, 0, 0.25])
		ax.set_xlabel(r'QSO rest-frame wavelength [$\AA$]', fontsize=18)
		ax.set_ylabel('Residuals', fontsize=18)
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=14)
		ax.yaxis.set_label_coords(-0.06, 0.5)
		# ax.legend(fontsize = 16)
		plt.subplots_adjust(hspace=0)
		ax = plt.gca()
		ax.tick_params(direction='in', which='both', labelsize=16)
		if save_plot == True:
			if directory is None or filename is None:
				if directory is None:
					print('Directory is not defined, final fit plot will not be saved.')
				elif filename is None:
					print('Filename not specified, final fit plot will not be saved.')
			else:
				plt.savefig(join(directory, filename), dpi=300, format='pdf',
				            transparent=True, bbox_inches='tight', pad_inches=0)
				print('Final fit plot succesfully saved.')
		if show_plot == True:
			plt.show()
		else:
			plt.close()

	def _bin_data(self):

		self.wave = np.array([np.mean(self.wave[n:n + self.binning]) for n in
		                      range(0, len(self.wave), self.binning)])
		self.flux = np.array([np.mean(self.flux[n:n + self.binning]) for n in
		                      range(0, len(self.flux), self.binning)])
		self.std = np.array(
			[np.mean(self.std[n:n + self.binning]) / np.sqrt(self.binning) for n in
			 range(0, len(self.std), self.binning)])
