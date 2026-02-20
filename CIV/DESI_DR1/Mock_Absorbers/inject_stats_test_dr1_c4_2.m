%% Load Data
load('/home/cosmic/IGM-Mining-James/output/2-17-26_mock_absorbers_injection/processed/processed_sigma_125_Spline_2200-2300_DESI-dr1.mat')
load('/home/cosmic/IGM-Mining-James/output/2-17-26_mock_absorbers_injection/processed/preloaded_qsos_DESI-dr1.mat')
load('/home/cosmic/IGM-Mining-James/output/2-17-26_mock_absorbers_injection/processed/SALSA_catalog.mat')

%% Setup tables of absorbers + data
% all mock absorber data are in 1000X5 arrays
mock_abs_z      = cat(2,sim_z_15,sim_z_20,sim_z_30,sim_z_40,sim_z_50); % simulated catalog z values for each redshift
mock_abs_EW1548 = cat(2,sim_EW1548_15,sim_EW1548_20,sim_EW1548_30,... % simulated rest EW values
    sim_EW1548_40,sim_EW1548_50);
mock_abs_EW1550 = cat(2,sim_EW1550_15,sim_EW1550_20,sim_EW1550_30,... % simulated rest EW values
    sim_EW1550_40,sim_EW1550_50);
mock_abs_N1548 = cat(2,sim_N1548_15,sim_N1548_20,sim_N1548_30,... % simulated column density values
    sim_N1548_40,sim_N1548_50);
mock_abs_N1550 = cat(2,sim_N1550_15,sim_N1550_20,sim_N1550_30,... % simulated column density values
    sim_N1550_40,sim_N1550_50);

% all detected absorber data are in 1000x7 arrays
det_abs_mask = (p_c4>0.85); %filter to only >85% prob absorbers
det_map_z = map_z_c4L2.*det_abs_mask;
det_map_z(det_map_z == 0) = NaN;
REW_1548_dr1 = REW_1548_dr1(1:1000,7); 
REW_1550_dr1 = REW_1550_dr1(1,1000:7);
det_map_REW1548 = REW_1548_dr1.*det_abs_mask;
det_map_REW1550 = REW_1550_dr1.*det_abs_mask;
det_map_N1548 = map_N_c4L2.*det_abs_mask;


% Filter mock absorbers to ones in range
valid_abs   = masked_abs_pixels == 0; % injected absorbers that do not have masked pixels
mock_abs_z          = mock_abs_z.* valid_abs; %index by valid mock absorbers
mock_abs_z(mock_abs_z == 0) = NaN;              %set 0s to NaNs
mock_abs_EW1548     = mock_abs_EW1548.* valid_abs;
mock_abs_EW1548(mock_abs_z == 0) = NaN;
mock_abs_REW1548 = mock_abs_EW1548./(1+mock_abs_z); %convert EW to REW
mock_abs_EW1550     = mock_abs_EW1550.* valid_abs;
mock_abs_EW1550(mock_abs_z == 0) = NaN;
mock_abs_REW1548 = mock_abs_EW1550./(1+mock_abs_z); %convert EW to REW
mock_abs_N1548      = mock_abs_N1548.* valid_abs;
mock_abs_N1548(mock_abs_z == 0) = NaN;
mock_abs_N1550      = mock_abs_N1550.* valid_abs;
mock_abs_N1550(mock_abs_z == 0) = NaN;

%% matching detected absorbers to mock absorbers
[num_quasars,num_c4] = size(p_c4);
det_abs_z = NaN(size(mock_abs_z));
det_abs_N1548 = NaN(size(mock_abs_z));
det_abs_SNR = NaN(size(mock_abs_z));
det_abs_REW1548 = NaN(size(mock_abs_z));
det_abs_REW1550 = NaN(size(mock_abs_z));
for i = 1:num_quasars
    for j = 1:5 %number of mock redshift bins
        this_mock_abs_z = mock_abs_z(i,j);
        if this_mock_abs_z==0   %skip if there is no valid mock absorber
            continue
        end
        %check for detected redshifts close to current mock redshift
        abs_dif = abs(det_map_z(i,:)-this_mock_abs_z);
        [~,jj] = find(abs_dif<0.1);
        if isempty(jj)
            continue
        end
        det_abs_z(i,j)      = det_map_z(i,jj);
        det_abs_N1548(i,j)  = det_map_N1548(i,jj);
        det_abs_SNR(i,j)    = det_map_SNR(i);
        % det_abs_EW1548 = det_map_EW1548[i,jj];
        % det_abs_EW1548 = det_map_EW1548[i,jj];
    end
end

%% Reshape arrays to vectors
% 1000X5 arrays --> 5000X1 vectors
mock_abs_z          = reshape(mock_abs_z,[],1);
mock_abs_REW1548     = reshape(mock_abs_REW1548,[],1);
mock_abs_REW1550     = reshape(mock_abs_REW1550,[],1);
mock_abs_N1548      = reshape(mock_abs_N1548,[],1);
mock_abs_N1550      = reshape(mock_abs_N1550,[],1);
det_abs_z           = reshape(det_abs_z,[],1);
det_abs_N1548       = reshape(det_abs_N1548,[],1);
det_abs_SNR         = reshape(det_abs_SNR,[],1);

%clear empty rows
mock_abs_bin = ~isnan(mock_abs_z);
mock_abs_z          = mock_abs_z(mock_abs_bin);
mock_abs_REW1548     = mock_abs_REW1548(mock_abs_bin);
mock_abs_REW1550     = mock_abs_REW1550(mock_abs_bin);
mock_abs_N1548      = mock_abs_N1548(mock_abs_bin);
mock_abs_N1550      = mock_abs_N1550(mock_abs_bin);
det_abs_z           = det_abs_z(mock_abs_bin);
det_abs_N1548       = det_abs_N1548(mock_abs_bin);
det_abs_SNR         = det_abs_SNR(mock_abs_bin);

det_abs_bin = ~isnan(det_abs_z); %create binary array of where our detected absorbers are

%% make completeness curve
figure(1)
REWrange = (12:0.5:17);
for SN = 2:8
    SNR_ind = (det_abs_SNR >= SN);
    for REWi = 1:numel(REWrange)
        REW = REWrange(REWi);
        REW_ind = (mock_abs_REW1548 >= REW);
        this_TP = nnz((mock_abs_bin(SNR_ind)) & (det_abs_bin(SNR_ind)));