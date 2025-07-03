% read_spec: loads data from our DESI dr1 data fits files;

function [wavelengths, flux, noise_variance, pixel_mask, sigma_pixel] = read_spec_DESI(filename)

  % mask bits to consider
  BRIGHTSKY = 24;

  measurements = fitsread(filename, ...
          'binarytable',  1, ...
          'tablecolumns', 1:3);

  % coadded calibrated flux  10^-17 erg s^-1 cm^-2 A^-1
  wavelengths = measurements{1};
  flux = measurements{2};
  inverse_noise_variance = measurements{3};

  % sigma of pixel in units of number of pixel 
  sigma_pixel = ones(size(flux))*0.94147; %Estimate, test value

  % derive noise variance
  noise_variance = 1 ./ (inverse_noise_variance);

  % derive bad pixel mask, remove pixels considered very bad
  % (FULLREJECT, NOSKY, NODATA); additionally remove pixels with
  % BRIGHTSKY set
  pixel_mask = ...
      (inverse_noise_variance <= 0) |...
      (noise_variance <=0) | (sigma_pixel<=0) | isnan(sigma_pixel) | ...
      (isnan(noise_variance) | isnan(inverse_noise_variance) | isinf(sigma_pixel) | ...
      isinf(noise_variance) | isinf(inverse_noise_variance)); %| bitget(and_mask, BRIGHTSKY)

end