% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
tic;

% load(sprintf('%s/catalog', processed_directory(releaseTest)), ...
%     variables_to_load{:});

num_quasars = numel(all_zqso_dr1);

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_sigma_pixel     =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);


% to track reasons for filtering out QSOs
filter_flags = zeros(num_quasars, 1, 'uint8');

%BAL filter not needed, applied in data download

% coverage : Shouldn't be needed, also applied in data download
%ind = (1310*(1+all_zqso_dr16)<3650) | (1548*(1+all_zqso_dr16)>10000); 
%filter_flags(ind)  = bitset(filter_flags(ind), 2, true);

num_quasars=100;
for i = 1:num_quasars

  if (filter_flags(i)~=0)
    continue;
  end

  %------dr16----------
  
    [this_wavelengths, this_flux, this_noise_variance, this_pixel_mask, this_sigma_pixel] ...
    = file_loader_DESI(all_QSO_ID_dr1(i));
   % % Masking Sky lines 

  this_pixel_mask((abs(this_wavelengths-5579)<5) & (abs(this_wavelengths-6302)<5))=1;
  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, all_zqso_dr16(i));
  % normalize flux

  ind = (this_rest_wavelengths >= normalization_min_lambda) & ...
        (this_rest_wavelengths <= normalization_max_lambda) & ...
        (~this_pixel_mask);


  this_median = median(this_flux(ind), "omitmissing");

  % bit 2: cannot normalize (all normalizing pixels are masked)

  if (isnan(this_median))
    filter_flags(i) = bitset(filter_flags(i), 3, true);
    continue;
  end

  ind = (this_rest_wavelengths >= min_lambda) & ...
          (this_rest_wavelengths <= max_lambda); % & ...
          (~this_pixel_mask);

  % bit 3: not enough pixels available
  if (nnz(ind) < min_num_pixels)
    filter_flags(i) = bitset(filter_flags(i), 4, true);
    continue;
  end

  all_normalizers(i) = this_median;

  this_flux           = this_flux           / this_median;
  this_noise_variance = this_noise_variance / this_median^2;

  % add one pixel on either side
  available_ind = find(~ind & ~this_pixel_mask);
  ind(min(available_ind(available_ind > find(ind, 1, 'last' )))) = true;
  ind(max(available_ind(available_ind < find(ind, 1, 'first')))) = true;

  all_wavelengths{i}    =    this_wavelengths(ind);%no longer (ind)
  all_flux{i}           =           this_flux(ind);
  all_noise_variance{i} = this_noise_variance(ind);
  all_pixel_mask{i}     =     this_pixel_mask(ind);
  all_sigma_pixel{i}    =     this_sigma_pixel(ind); 


  fprintf('loaded quasar %i of %i (%i)\n', ...
          i, num_quasars, all_QSO_ID_dr1(i));

end

variables_to_save = {'loading_min_lambda', 'loading_max_lambda', ...
                     'normalization_min_lambda', 'normalization_max_lambda', ...
                     'min_num_pixels', 'all_wavelengths', 'all_flux', ...
                     'all_noise_variance', 'all_pixel_mask', ...
                     'all_normalizers', 'all_sigma_pixel', 'filter_flags'};
save(sprintf('%s/preloaded_qsos_%s', processed_directory(releaseTest), testing_set_name), ...
     variables_to_save{:}, '-v7.3');

% write new filter flags to catalog
save(sprintf('%s/filter_flags', processed_directory(releaseTest)), ...
     'filter_flags', '-v7.3');

toc;

