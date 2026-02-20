%% Load Data
clear
clc
load('/home/cosmic/IGM-Mining-James/output/2-20-26_mock_absorbers_injection/processed/processed_sigma_125_Spline_2200-2300_DESI-dr1.mat')
load('/home/cosmic/IGM-Mining-James/output/2-20-26_mock_absorbers_injection/processed/preloaded_qsos_DESI-dr1.mat')
load('/home/cosmic/IGM-Mining-James/output/2-20-26_mock_absorbers_injection/processed/SALSA_catalog.mat')

%% Setup
z_centers = [1.5 2.0 3.0 4.0 5.0]; 
nqso      = size(masked_abs_pixels,1);
nbin_z    = numel(z_centers);
maxcol    = 7;
EW_bins = [0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.80 0.95 1.25 1.55 2.0 3.0 4.0 5.0 inf];
N_bins  = [11:18,inf];
nbin_EW = numel(EW_bins)-1;
nbin_N  = numel(N_bins)-1;

valid_abs   = masked_abs_pixels == 0; % injected absorbers that do not have masked pixels
detected    = p_c4 > 0.85; % absorbers that we detect
detected_zc4  = map_z_c4L2 .* detected; % z values of detected absorber
detected_N  = map_N_c4L2 .* detected;
% REW_1548_dr1 = REW_1548_dr1;
% REW_1550_dr1 = REW_1550_dr1;
detected_REW1548 = REW_1548_dr1 .* detected;
detected_REW1550 = REW_1550_dr1 .* detected;

sim_z      = cat(2,sim_z_15,sim_z_20,sim_z_30,sim_z_40,sim_z_50); % simulated catalog z values for each redshift
sim_REW1548 = cat(2,sim_EW1548_15,sim_EW1548_20,sim_EW1548_30,... % simulated rest EW values
                   sim_EW1548_40,sim_EW1548_50)./(1 +sim_z);
sim_REW1550 = cat(2,sim_EW1550_15,sim_EW1550_20,sim_EW1550_30,... % simulated rest EW values
                   sim_EW1550_40,sim_EW1550_50)./(1 +sim_z);
sim_N1548 = cat(2,sim_N1548_15,sim_N1548_20,sim_N1548_30,... % simulated column density values
    sim_N1548_40,sim_N1548_50);
colors = ['r','g','b','c','m','k'];
SNR_range = 3:8; % SNR < 3 are already cut out in processes
z_det_all = []; z_cat_all = []; N_det_all = []; N_cat_all = [];
REW1548_det_all = []; REW1550_det_all = []; REW1548_cat_all = []; REW1550_cat_all = [];
QSO_ind = []; all_completeness_REW = []; all_completeness_N = [];

for j = 1:numel(SNR_range)
    SNR_cut    = all_SNR > SNR_range(j); % SNR cut to look at completness for higher SNR cuts
    
    
    
    % Z catalog vs Z detected
    detected_abs_logical = false(nqso,nbin_z); % initialize logical array
    
    for iz = 1:nbin_z % iterates through all redshift columns from the SALSA catalog for c4
        in_z      = abs(detected_zc4 - z_centers(iz)) < 0.1; % find detected abs for current redshift
        mask      = in_z(:,1:maxcol) & valid_abs(:,iz); % apply row of valid_abs for current redshift
        % to mask and detected abs that are not valid logical array to all rows of in_z

        [row,~]   = find(mask); % rows in both sim and detected are the same spectra, 
        % finds rows with found absorber for current redshift (should only be 1)

        detected_abs_logical(row,iz) = true; % logical array used to mask simulated absorbers
        % that did not get found, created one row at a time for each current redshift

        if j == 1 % only want this for all data (before we start cutting low SNR spectra)
            z_det_all = [z_det_all; detected_zc4(mask)]; % all redshift values our code finds for detected abs
            z_cat_all = [z_cat_all; sim_z(row,iz)]; % all redshift values from SALSA for detected abs
            N_det_all = [N_det_all; detected_N(mask)];
            N_cat_all = [N_cat_all; sim_N1548(row,iz)];
            REW1548_det_all = [REW1548_det_all; detected_REW1548(mask)];
            REW1548_cat_all = [REW1548_cat_all; sim_REW1548(row,iz)];
            REW1550_det_all = [REW1550_det_all; detected_REW1550(mask)];
            REW1550_cat_all = [REW1550_cat_all; sim_REW1550(row,iz)];
            QSO_ind = [QSO_ind;row];
        end
        
    end
    
    
    % EW completeness
    detected_abs_logical = detected_abs_logical & SNR_cut; % mask absorbers below SNR cut for current loop
    all_EW = sim_REW1548(1:1850,:) .* valid_abs .* SNR_cut; % catalog absobers with only valid and SNR masks
    all_N = sim_N1548(1:1850,:) .* valid_abs .* SNR_cut; % catalog absobers with only valid and SNR masks
    
    completeness_REW = nan(1,nbin_EW);
    completeness_N = nan(1,nbin_N);

    EW_all_vec   = all_EW(:); % 1000 x 5 -> 5000 x 1 vector
    det_all_vec  = detected_abs_logical(:);
    N_all_vec = all_N(:);
    
    for i = 1:nbin_EW % completness for each EW bin
        ind = EW_all_vec >= EW_bins(i) & EW_all_vec < EW_bins(i+1);
        completeness_REW(i) = nnz(det_all_vec(ind)) / nnz(ind);
    end
    all_completeness_REW = [all_completeness_REW; completeness_REW];
    % for i = 1:nbin_N % completness for each N bin
    %     ind_N = log10(N_all_vec) >= N_bins(i) & log10(N_all_vec) < N_bins(i+1);
    %     completeness_N(i) = nnz(det_all_vec(ind_N)) / nnz(ind_N);
    % end
    % all_completeness_N = [all_completeness_N; completeness_N];

end

%% Plots
% figure(1) % Column Density Completness
% for j = 1:numel(SNR_range)
%     plot(1:nbin_N,all_completeness_N(j,:),'-o','Color',colors(j),...
%          'DisplayName', sprintf('SNR > %.2f', SNR_range(j)))
%     hold on
% end
% hold off
% xlabel('$\mathrm{N}_{\mathrm{C\,IV},1548}\ [\AA]$','Interpreter','latex')
% ylabel('completeness')
% xticks(1:nbin_N)
% xticklabels({'11-12','12-13','13-14','14-15','15-16','16-17','17-18','18+'})


figure(2) % REW Completness
for j = 1:numel(SNR_range)
    plot(1:nbin_EW,all_completeness_REW(j,:),'-o','Color',colors(j),...
         'DisplayName', sprintf('SNR > %.2f', SNR_range(j)))
    hold on
end
hold off
xlabel('$\mathrm{EW}^{\mathrm{rest}}_{\mathrm{C\,IV},1548}\ [\AA]$','Interpreter','latex')
ylabel('completeness')
xticks(1:nbin_EW)
xticklabels({'0.05-0.15','0.15-0.25','0.25-0.35','0.35-0.45','0.45-0.55',...
    '0.55-0.65','0.65-0.8','0.8-0.95','0.95-1.25','1.25-1.55','1.55-2.0',...
    '2.0-3.0','3.0-4.0','4.0-5.0','>5.0'})
xtickangle(90)
legend(Location="northwest")

fitfun = @(a,b,x) a*x + b; % linear fit function
z_fit  = fit(z_det_all,z_cat_all,fitfun);
%%
figure(3) % z_detected vs. z_catalog
plot(z_det_all,z_cat_all,'.r'); 
hold on
x_array = linspace(0,2.2,1000);
plot(x_array,fitfun(z_fit.a,z_fit.b,x_array),'k--')
xlabel('z_{detected}'); ylabel('z_{catalog}')
hold off

% N_fit = fit(log10(N_cat_all),N_det_all,fitfun);
% ci_N = confint(N_fit,0.68);
x2_array = linspace(12,19,1000);

figure(4) % N_catalog vs. N_detected
scatter(log10(N_cat_all),N_det_all,'b.');
hold on
% plot(x2_array,fitfun(N_fit.a,N_fit.b,x2_array),'k--')
% plot(x2_array,fitfun(ci_N(1,1),N_fit.b,x2_array),'k:')
% plot(x2_array,fitfun(ci_N(2,1),N_fit.b,x2_array),'k:')
plot (x2_array,x2_array,'r--')
hold off
ylabel('log(N_{detected})')
xlabel('log(N_{catalog})')
axis equal

% REW_fit = fit(REW1548_cat_all,REW1548_det_all,fitfun);
% ci_REW = confint(REW_fit,0.68);
x4_array = linspace(0,4,100);

figure(5) % REW catalog vs REW detected
plot(REW1548_cat_all,REW1548_det_all,'.r'); 
hold on
% plot(x4_array,fitfun(REW_fit.a,REW_fit.b,x4_array),'r--')
% plot(x4_array,fitfun(ci_REW(1,1),REW_fit.b,x4_array),'r:')
% plot(x4_array,fitfun(ci_REW(2,1),REW_fit.b,x4_array),'r:')
plot(0:4,0:4,'b--')
xlabel('REW_{cat}'); ylabel('REW_{det}')
% legend('data',sprintf('y = (%.3f%c%.3f)x + %.3f',REW_fit.a,...
%     177,abs(REW_fit.a-ci_REW(1,1)),REW_fit.b))
hold off

figure(6) % Doublet Ratio detected vs. DR cat
plot(REW1548_det_all./REW1550_det_all,REW1548_cat_all./REW1550_cat_all,'.r'); 
hold on
plot(1:2,1:2,'b--')
xlabel('DR_{detected}'); ylabel('DR_{catalog}')
hold off


%% More Plots

figure(7) % Delta REW for all QSO ind
scatter(QSO_ind,REW1548_cat_all-REW1548_det_all,'k.')
xlabel('QS0 Index')
ylabel('\DeltaEW')

% delta_fit = fit(REW1548_cat_all,REW1548_cat_all-REW1548_det_all,fitfun);
% x3_array = linspace(0,3,1000);

figure(8) % Delta REW as a function of REW_catalog
plot(REW1548_cat_all,REW1548_cat_all-REW1548_det_all,'k.')
hold on
% plot(x3_array, fitfun(delta_fit.a,delta_fit.b,x3_array),'b--')
hold off
xlabel('REW_{catalog}')
ylabel('REW_{catalog} - REW_{detected}')
% legend('data',sprintf('y = %.3fx + %.3f',delta_fit.a,delta_fit.b))

figure(9) % Distribution of sig and N samples
scatter(map_sigma_c4L2(:),map_N_c4L2(:))
xlabel('\sigma')
ylabel('log(N)')