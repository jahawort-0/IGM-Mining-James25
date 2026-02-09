% Build catalogs usable for spectra from DESI dr1 and from SALSA simulated
% absorbers

DESI_dr1q = ...
fitsread('/home/cosmic/desi-env/sim_abs_CIV/QSO_catalog.fits', 'binarytable');
all_tid_dr1               = DESI_dr1q{1};
all_tid_dr1               = all_tid_dr1(null_search_filter);
all_tid_dr1               = all_tid_dr1(null_QSO_ind);
all_zqso_dr1              = DESI_dr1q{2};
all_zqso_dr1              = all_zqso_dr1(null_search_filter);
all_zqso_dr1              = all_zqso_dr1(null_QSO_ind);
all_RA_dr1                = DESI_dr1q{3};
all_RA_dr1                = all_RA_dr1(null_search_filter);
all_RA_dr1                = all_RA_dr1(null_QSO_ind);
all_DEC_dr1               = DESI_dr1q{4};
all_DEC_dr1               = all_DEC_dr1(null_search_filter);
all_DEC_dr1               = all_DEC_dr1(null_QSO_ind);

% Load all 5 simulated absorber catalogs 
sim_cat_15 = ... % z = 1.5
fitsread('/home/cosmic/desi-env/sim_abs_CIV/CIV_Simulated_Catalog_z1.5.fits', 'binarytable');
sim_id_15                 = sim_cat_15{1};
sim_z_15                  = sim_cat_15{2};
sim_EW1548_15             = sim_cat_15{3};
sim_EW1550_15             = sim_cat_15{4};
sim_N1548_15              = sim_cat_15{5};
sim_N1550_15              = sim_cat_15{6};

sim_cat_20 = ... % z = 2.0
fitsread('/home/cosmic/desi-env/sim_abs_CIV/CIV_Simulated_Catalog_z1.5.fits', 'binarytable');
sim_id_20                 = sim_cat_20{1};
sim_z_20                  = sim_cat_20{2};
sim_EW1548_20             = sim_cat_20{3};
sim_EW1550_20             = sim_cat_20{4};
sim_N1548_20              = sim_cat_20{5};
sim_N1550_20              = sim_cat_20{6};

sim_cat_30 = ... % z = 13.0
fitsread('/home/cosmic/desi-env/sim_abs_CIV/CIV_Simulated_Catalog_z1.5.fits', 'binarytable');
sim_id_30                 = sim_cat_30{1};
sim_z_30                  = sim_cat_30{2};
sim_EW1548_30             = sim_cat_30{3};
sim_EW1550_30             = sim_cat_30{4};
sim_N1548_30              = sim_cat_30{5};
sim_N1550_30              = sim_cat_30{6};

sim_cat_40 = ... % z = 4.0
fitsread('/home/cosmic/desi-env/sim_abs_CIV/CIV_Simulated_Catalog_z1.5.fits', 'binarytable');
sim_id_40                 = sim_cat_40{1};
sim_z_40                  = sim_cat_40{2};
sim_EW1548_40             = sim_cat_40{3};
sim_EW1550_40             = sim_cat_40{4};
sim_N1548_40              = sim_cat_40{5};
sim_N1550_40              = sim_cat_40{6};

sim_cat_50 = ... % z = 5.0
fitsread('/home/cosmic/desi-env/sim_abs_CIV/CIV_Simulated_Catalog_z1.5.fits', 'binarytable');
sim_id_50                 = sim_cat_50{1};
sim_z_50                  = sim_cat_50{2};
sim_EW1548_50             = sim_cat_50{3};
sim_EW1550_50             = sim_cat_50{4};
sim_N1548_50              = sim_cat_50{5};
sim_N1550_50              = sim_cat_50{6};

% save catalog 
variables_to_save = {'all_tid_dr1', 'all_RA_dr1', 'all_DEC_dr1', 'all_zqso_dr1'};
save(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3');

variables_to_save = {'sim_id_15','sim_z_15','sim_EW1548_15','sim_EW1550_15', 'sim_N1548_15','sim_N1550_15'...
    'sim_id_20','sim_z_20','sim_EW1548_20','sim_EW1550_20', 'sim_N1548_20','sim_N1550_20'...
    'sim_id_30','sim_z_30','sim_EW1548_30','sim_EW1550_30', 'sim_N1548_30','sim_N1550_30'...
    'sim_id_40','sim_z_40','sim_EW1548_40','sim_EW1550_40', 'sim_N1548_40','sim_N1550_40'...
    'sim_id_50','sim_z_50','sim_EW1548_50','sim_EW1550_50', 'sim_N1548_50','sim_N1550_50'};
save(sprintf('%s/SALSA_catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3')