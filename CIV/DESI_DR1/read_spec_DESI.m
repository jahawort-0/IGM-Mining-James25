% read_spec: loads data from our DESI dr1 data fits files;

function [wavelengths, flux, noise_variance, pixel_mask, sigma_pixel] = read_spec_DESI(filename)

  % mask bits to consider
  %BRIGHTSKY = 24;

  measurements = fitsread(filename, ...
          'binarytable',  1, ...
          'tablecolumns', 1:5);

  %Data values from .fits files
  wavelengths = measurements{1};
  flux = measurements{2};
  inverse_noise_variance = measurements{3};
  default_pixel_mask = measurements{4};
  sigma_pixel = measurements{5};

  % derive noise variance
  noise_variance = 1 ./ (inverse_noise_variance);

  % find low signal to noise pixels (may be cosmic rays etc that have been
  % missed)
  signalnoise = flux./noise_variance;
  snflag = signalnoise<2;

  % derive bad pixel mask, remove pixels considered very bad
  % (FULLREJECT, NOSKY, NODATA); additionally remove pixels with
  % BRIGHTSKY set
  pixel_mask = ...
      (inverse_noise_variance <= 0) |(noise_variance <=0) | (sigma_pixel<=0) |...
      isnan(sigma_pixel) | (isnan(noise_variance) | isnan(inverse_noise_variance) |...
      isinf(sigma_pixel) | isinf(noise_variance) | isinf(inverse_noise_variance)|...
      default_pixel_mask)| snflag; %| bitget(and_mask, BRIGHTSKY)

end