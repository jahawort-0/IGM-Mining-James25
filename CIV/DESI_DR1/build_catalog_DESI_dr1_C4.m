% Build catalogs usable for spectra from dr16

DESI_dr1q = ...
fitsread('QSO_catalog.fits', 'binarytable');
all_tid_dr1               = DESI_dr1q{1};
all_zqso_dr1              = DESI_dr1q{2};
all_RA_dr1                = DESI_dr1q{3};
all_DEC_dr1               = DESI_dr1q{4};


num_quasars_dr1 = numel(all_zqso_dr1);%= 1000;
all_QSO_ID_dr1  = all_tid_dr1;


% save catalog 
variables_to_save = {'all_QSO_ID_dr1', 'all_RA_dr1', 'all_DEC_dr1', 'all_zqso_dr1'};
save(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3');