% set_parameters: sets various parameters for the CIV detection 
% pipeline
% Designed for using DR16 spectra  
%flags for changes
extrapolate_subdla = 0; %0 = off, 1 = on
add_proximity_zone = 0;
integrate          = 1;
optTag = [num2str(integrate), num2str(extrapolate_subdla), num2str(add_proximity_zone)];

% physical constants
civ_1548_wavelength = 1548.1949462890625;		 % CIV transition wavelength  Å
civ_1550_wavelength =  1550.77001953125; 		 % CIV transition wavelength  Å
speed_of_light = 299792458;                   % speed of light                     m s⁻¹

% converts relative velocity in km s^-1 to redshift difference
kms_to_z = @(kms) (kms * 1000) / speed_of_light;

% utility functions for redshifting
emitted_wavelengths = ...
    @(observed_wavelengths, z) (observed_wavelengths / (1 + z));

observed_wavelengths = ...
    @(emitted_wavelengths,  z) ( emitted_wavelengths * (1 + z));

% %-------dr7
% release = 'dr7';
% % download Cooksey's dr7 spectra from this page: 
% % http://www.guavanator.uhh.hawaii.edu/~kcooksey/SDSS/CIV/index.html 
% % go to table: "SDSS spectra of the sightlines surveyed for C IV."
% file_loader = @(mjd, plate, fiber_id) ...
%   (read_spec_dr7(sprintf('data/dr7/spectro/1d_26/%04i/1d/spSpec-%05i-%04i-%03i.fit',...
%   plate, mjd, plate, fiber_id)));
releaseTest='dr1';
releasePrior='dr7'; %Need to decide on a prior dataset


% file loading parameters
loading_min_lambda = 1310;          % range of rest wavelengths to load  Å
loading_max_lambda = 1555;                    
% The maximum allowed is set so that even if the peak is redshifted off the end, the
% quasar still has data in the range

% preprocessing parameters
%z_qso_cut      = 2.15;                   % filter out QSOs with z less than this threshold
%z_qso_cut      = 1.7;                    % according to Cooksey z>1.7                      
min_num_pixels = 600;                     % minimum number of non-masked pixels NEED ADJUSTMENT

% normalization parameters
normalization_min_lambda = 1420;                    % range of rest wavelengths to use for flux norm Å
normalization_max_lambda = 1475; 

% null model parameters             XXX DESI has a 0.8A wavelength spacing
min_lambda         = loading_min_lambda+1;                   % range of rest wavelengths to       Å
max_lambda         = loading_max_lambda-1;                    %   model
dlambda            = 0.8;                    % separation of wavelength grid      Å
k                  = 20;                      % rank of non-diagonal contribution
max_noise_variance = 0.5^2;                   % maximum pixel noise allowed during model training
h                  = 2;                     % masking par to remove CIV region 
nAVG               = 20;                     % number of points added between two 
                                            % observed wavelengths to make the Voigt finer
% optimization parameters
minFunc_options =               ...           % optimization options for model fitting
    struct('MaxIter',     10000, ...
           'MaxFunEvals', 10000);

% C4 model parameters: parameter samples (for Quasi-Monte Carlo)
num_C4_samples           = 50000;                  % number of parameter samples
alpha                    = 0.9;                    % weight of KDE component in mixture
uniform_min_log_nciv     = 12.5;                   % range of column density samples    [cm⁻²]
uniform_max_log_nciv     = 16.1;                   % from uniform distribution
fit_min_log_nciv         = uniform_min_log_nciv;                   % range of column density samples    [cm⁻²]
fit_max_log_nciv         = uniform_max_log_nciv;                   % from fit to log PDF

min_sigma                = 5e5;                   % cm/s -> b/sqrt(2) -> min Doppler par from Cooksey
max_sigma                = 115e5;                   % cm/s -> b/sqrt(2) -> max Doppler par from Cooksey
vCut                     = 3000;                    % maximum cut velocity for CIV system 
RejectionSampling        = 0;
% model prior parameters
prior_z_qso_increase = kms_to_z(30000);       % use QSOs with z < (z_QSO + x) for prior

% instrumental broadening parameters
width = 3;                                    % width of Gaussian broadening (# pixels)
pixel_spacing = 1e-4;                         % wavelength spacing of pixels in dex

% DLA model parameters: absorber range and model
num_lines = 2;                                % number of members of CIV series to use

max_z_cut = kms_to_z(vCut);                   % max z_DLA = z_QSO - max_z_cut
% max_z_c4 = @(wavelengths, z_qso) ...         % determines maximum z_DLA to search
%     (max(wavelengths)/civ_1548_wavelength - 1) - max_z_cut;
max_z_c4 = @(z_qso, max_z_cut) ...         % determines maximum z_DLA to search
     z_qso - max_z_cut*(1+z_qso);
min_z_cut = kms_to_z(vCut);                   % min z_DLA = z_Ly∞ + min_z_cut
min_z_c4 = @(wavelengths, z_qso) ...         % determines minimum z_DLA to search
    max(min(wavelengths) / civ_1548_wavelength - 1,                          ...
        observed_wavelengths(1310, z_qso) / civ_1548_wavelength - 1 + ...
        min_z_cut);
% min_z_c4 = @(wavelengths, z_qso) ...         % determines minimum z_DLA to search
%      min(wavelengths) / civ_1548_wavelength - 1;
train_ratio =1;
training_set_name = 'C13-full';  %% Change
testing_set_name = 'DESI-dr1';  
max_civ = 7; 
% base directory for all data 
%Change to something new on Cosmic
base_directory = 'output';
% utility functions for identifying various directories
distfiles_directory = @(release) ...
   sprintf('%s/%s/distfiles', base_directory, release);

spectra_directory   = @(release)...
   sprintf('%s/%s/sas/dr1/sdss/spectro/redux/v5_13_0/spectra/lite', base_directory, release);

processed_directory = @(release) ...
   sprintf('%s/%s/processed', base_directory, release);

c4_catalog_directory = @(name) ...
   sprintf('%s/C4_catalogs/%s/processed', base_directory, name);

   
% replace with @(varargin) (fprintf(varargin{:})) to show debug statements
% fprintf_debug = @(varargin) (fprintf(varargin{:}));
% fprintf_debug = @(varargin) ([]);

%---------DESI dr1-------------

file_loader_DESI = @(tid)...
    (read_spec_DESI(sprintf("/home/cosmic/DESI_data/%s.fits",tid)));
    