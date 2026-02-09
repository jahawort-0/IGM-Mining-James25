% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
tic;

load(sprintf('%s/catalog', processed_directory(releaseTest)), ...
   variables_to_load{:});
load(sprintf('%s/SALSA_catalog', processed_directory(releaseTest)), ...
   variables_to_load{:});

% load all optical depth data
data = fitsread('/home/cosmic/desi-env/sim_abs_CIV/Tau_CIV_z1.5.fits');
wave = data(1,:);
tau_data_15 = data(2:end,:);
data = fitsread('/home/cosmic/desi-env/sim_abs_CIV/Tau_CIV_z2.0.fits');
tau_data_20 = data(2:end,:);
data = fitsread('/home/cosmic/desi-env/sim_abs_CIV/Tau_CIV_z3.0.fits');
tau_data_30 = data(2:end,:);
data = fitsread('/home/cosmic/desi-env/sim_abs_CIV/Tau_CIV_z4.0.fits');
tau_data_40 = data(2:end,:);
data = fitsread('/home/cosmic/desi-env/sim_abs_CIV/Tau_CIV_z5.0.fits');
tau_data_50 = data(2:end,:);

%num_quasars = numel(all_zqso_dr1);     %number also comes from run_script

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_sigma_pixel    =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);
all_optical_depth  =  cell(num_quasars, 1);
all_isAbs          =  cell(num_quasars, 5);



% to track reasons for filtering out QSOs
filter_flags = zeros(num_quasars, 1, 'uint8');

%BAL filter not needed, applied in data download

% coverage : Shouldn't be needed, also applied in data download
%ind = (1310*(1+all_zqso_dr16)<3650) | (1548*(1+all_zqso_dr16)>10000); 
%filter_flags(ind)  = bitset(filter_flags(ind), 2, true);

%num_quasars=29;
for i = 1:num_quasars 

  if (filter_flags(i)~=0)
    continue;
  end

  %------dr16----------
  
    [this_wavelengths, this_flux, this_noise_variance, this_pixel_mask, this_sigma_pixel] ...
    = file_loader_DESI(all_tid_dr1{i});
    % Interpolate the tau data to the wavelegnth space from DESI
    this_tau_interp_15 = interp1(wave', tau_data_15(i,:)',this_wavelengths, 'pchip');  
    this_tau_interp_20 = interp1(wave', tau_data_20(i,:)',this_wavelengths, 'pchip');  
    this_tau_interp_30 = interp1(wave', tau_data_30(i,:)',this_wavelengths, 'pchip');  
    this_tau_interp_40 = interp1(wave', tau_data_40(i,:)',this_wavelengths, 'pchip');  
    this_tau_interp_50 = interp1(wave', tau_data_50(i,:)',this_wavelengths, 'pchip');  
    this_tau = exp(-this_tau_interp_15).*exp(-this_tau_interp_20)...
      .*exp(-this_tau_interp_30).*exp(-this_tau_interp_40).*exp(-this_tau_interp_50);
    this_flux = this_flux.*this_tau; % Convolved flux arrays for absorbers, with all 5 mock absorbers


    % Masking Sky lines 
    this_pixel_mask((abs(this_wavelengths-5579)<5) | (abs(this_wavelengths-6302)<5))=1;
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, all_zqso_dr1(i));
    % -- Normalization 
    ind = (this_rest_wavelengths >= normalization_min_lambda_3) & ...
           (this_rest_wavelengths <= normalization_max_lambda_3) & ...
           (~this_pixel_mask);
 
  this_median = median(this_flux(ind), "omitmissing");
  
  % bit 2: cannot normalize (all normalizing pixels are masked)
  if (isnan(this_median))
    filter_flags(i) = bitset(filter_flags(i), 3, true);
    continue;
  end

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

  all_wavelengths{i}    =  this_wavelengths;%no longer (ind)
  all_flux{i}           =  this_flux;
  all_noise_variance{i} =  this_noise_variance;
  all_pixel_mask{i}     =  this_pixel_mask;
  all_sigma_pixel{i}    =  this_sigma_pixel; 
  all_optical_depth{i}  =  this_tau;
  all_isAbs{i,1}        =  this_tau_interp_15 > 0;
  all_isAbs{i,2}        =  this_tau_interp_20 > 0;
  all_isAbs{i,3}        =  this_tau_interp_30 > 0;
  all_isAbs{i,4}        =  this_tau_interp_40 > 0;
  all_isAbs{i,5}        =  this_tau_interp_50 > 0;



  fprintf('loaded quasar %i of %i (%s)\n', ...
      i, num_quasars, all_tid_dr1{i});

end

variables_to_save = {'loading_min_lambda',...
                     'loading_max_lambda', ...
                     'min_num_pixels',...
                     'all_wavelengths',...
                     'all_flux', ...
                     'all_noise_variance',...
                     'all_pixel_mask', ...
                     'all_sigma_pixel',...
                     'filter_flags',...
                     'all_optical_depth',...
                     'all_isAbs'};
save(sprintf('%s/preloaded_qsos_%s.mat', processed_directory(releaseTest), testing_set_name), ...
     variables_to_save{:}, '-v7.3');

% write new filter flags to catalog
save(sprintf('%s/filter_flags', processed_directory(releaseTest)), ...
     'filter_flags', '-v7.3');

toc;