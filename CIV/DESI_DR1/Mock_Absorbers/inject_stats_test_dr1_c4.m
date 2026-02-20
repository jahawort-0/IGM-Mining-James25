%% Load Data
load('/home/cosmic/IGM-Mining-James/output/2-10-26_mock_absorbers_injection/processed/processed_sigma_125_Spline_2200-2300_DESI-dr1.mat')
load('/home/cosmic/IGM-Mining-James/output/2-10-26_mock_absorbers_injection/processed/preloaded_qsos_DESI-dr1.mat')
load('/home/cosmic/IGM-Mining-James/output/2-10-26_mock_absorbers_injection/processed/SALSA_catalog.mat')

%% Setup
z_centers = [1.5 2.0 3.0 4.0 5.0]; 
nqso      = size(masked_abs_pixels,1);
nbin_z    = numel(z_centers);
maxcol    = 7;
EW_bins = [0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.80 0.95 1.25 1.55 2.0 3.0 4.0 5.0 inf];
N_bins  = 11:18;
nbin_EW = numel(EW_bins)-1;
nbin_N  = numel(N_bins)-1;

valid_abs   = masked_abs_pixels == 0; % injected absorbers that do not have masked pixels
detected    = p_c4 > 0.85; % absorbers that we detect
detected_z  = map_z_c4L2 .* detected; % z values of detected absorbers
detected_N  = map_N_c4L2 .* detected;
% detected_EW1548 = map_REW_1548 .* detected;
% detected_EW1550 = map_REW_1550 .* detected;

sim_z      = cat(2,sim_z_15,sim_z_20,sim_z_30,sim_z_40,sim_z_50); % simulated catalog z values for each redshift
sim_EW1548 = cat(2,sim_EW1548_15,sim_EW1548_20,sim_EW1548_30,... % simulated rest EW values
                   sim_EW1548_40,sim_EW1548_50);
sim_EW2803 = cat(2,sim_EW1550_15,sim_EW1550_20,sim_EW1550_30,... % simulated rest EW values
                   sim_EW1550_40,sim_EW1550_50);
sim_N1548 = cat(2,sim_N1548_15,sim_N1548_20,sim_N1548_30,... % simulated column density values
    sim_N1548_40,sim_N1548_50);
sim_N1550 = cat(2,sim_N1550_15,sim_N1550_20,sim_N1550_30,... % simulated column density values
    sim_N1550_40,sim_N1550_50);
colors = ['r','g','b','c','m','k','y'];
figure(1)
SNR_range = 2:8; % SNR < 3 are already cut out in processes
z_det_all = []; z_cat_all = []; N_det_all = []; N_cat_all = [];
REW1548_det_all = []; REW1550_det_all = []; REW1548_cat_all = []; REW1550_cat_all = [];
for j = 1:numel(SNR_range)
    SNR_cut    = all_SNR > SNR_range(j); % SNR cut to look at completness for higher SNR cuts
    
    
    
    % Z catalog vs Z detected
    detected_abs_logical = false(nqso,nbin_z); % initialize logical array
    
    for iz = 1:nbin_z % iterates through all redshift columns from the SALSA catalog for MgII
        in_z      = abs(detected_z - z_centers(iz)) < 0.1; % find detected abs for current redshift
        mask      = in_z(:,1:maxcol) & valid_abs(:,iz); % apply row of valid_abs for current redshift
        % to mask and detected abs that are not valid logical array to all rows of in_z

        [row,~]   = find(mask); % rows in both sim and detected are the same spectra, 
        % finds rows with found absorber for current redshift (should only be 1)

        detected_abs_logical(row,iz) = true; % logical array used to mask simulated absorbers
        % that did not get found, created one row at a time for each current redshift

        if j == 1 % only want this for all data (before we start cutting low SNR spectra)
            z_det_all = [z_det_all; detected_z(mask)]; % all redshift values our code finds for detected abs
            z_cat_all = [z_cat_all; sim_z(row,iz)]; % all redshift values from SALSA for detected abs
            N_det_all = [N_det_all; detected_N(mask)];
            N_cat_all = [N_cat_all; sim_N1548(row,iz)];
            % REW1548_det_all = [REW1548_det_all; detected_EW1548(mask)];
            % REW1548_cat_all = [REW1548_cat_all; sim_EW1548(row,iz)];
            % REW1550_det_all = [REW1550_det_all; detected_EW1550(mask)];
            % REW1550_cat_all = [REW1550_cat_all; sim_EW1550(row,iz)];        
        end
        
    end
    

    
    % EW completeness
    detected_abs_logical = detected_abs_logical & SNR_cut; % mask absorbers below SNR cut for current loop
    % all_EW = sim_EW1548 .* valid_abs .* SNR_cut; % catalog absobers with only valid and SNR masks
    all_N = sim_N1548 .* valid_abs .* SNR_cut; % catalog absobers with only valid and SNR masks
    
    completeness = nan(1,nbin_EW);
    completeness_N = nan(1,nbin_N);

    % EW_all_vec   = all_EW(:); % 1000 x 5 -> 5000 x 1 vector
    det_all_vec  = detected_abs_logical(:);
    N_all_vec = all_N(:);
    
    % for i = 1:nbin_EW % completness for each EW bin
    %     ind = EW_all_vec >= EW_bins(i) & EW_all_vec < EW_bins(i+1);
    %     completeness(i) = nnz(det_all_vec(ind)) / nnz(ind);
    % end

    for i = 1:nbin_N % completness for each N bin
        ind_N = log10(N_all_vec) >= N_bins(i) & log10(N_all_vec) < N_bins(i+1);
        completeness_N(i) = nnz(det_all_vec(ind_N)) / nnz(ind_N);
    end
    
    plot(1:nbin_N,completeness_N,'-o','Color',colors(j),...
         'DisplayName', sprintf('SNR > %.2f', SNR_range(j)))
    hold on
end
hold off

%%
% xlabel('$\mathrm{EW}^{\mathrm{rest}}_{\mathrm{Mg\,II},2796}\ [\AA]$','Interpreter','latex')
xlabel('$\mathrm{N}_{\mathrm{C\,IC},1548}\ [\AA]$','Interpreter','latex')
ylabel('completeness')
% xticks(1:nbin_EW)
% xticklabels({'0.05-0.15','0.15-0.25','0.25-0.35','0.35-0.45','0.45-0.55',...
%     '0.55-0.65','0.65-0.8','0.8-0.95','0.95-1.25','1.25-1.55','1.55-2.0',...
%     '2.0-3.0','3.0-4.0','4.0-5.0','>5.0'})
xticks(1:nbin_N)
xticklabels({'11-12','12-13','13-14','14-15','15-16','16-17','17,-18'})
xtickangle(90)

legend(Location="northwest")

fitfun = @(a,b,x) a*x + b;
z_fit  = fit(z_det_all,z_cat_all,fitfun);

figure(2)
plot(z_det_all,z_cat_all,'.r'); 
hold on
x_array = linspace(0,2.2,1000);
plot(x_array,fitfun(z_fit.a,z_fit.b,x_array),'k--')
xlabel('z_{detected}'); ylabel('z_{catalog}')
hold off

N_fit = fit(log10(N_cat_all),N_det_all,fitfun);
ci_N = confint(N_fit,0.68);
x2_array = linspace(12,19,1000);
figure(3)
scatter(log10(N_cat_all),N_det_all,'b.');
hold on
plot(x2_array,fitfun(N_fit.a,N_fit.b,x2_array),'k--')
plot(x2_array,fitfun(ci_N(1,1),N_fit.b,x2_array),'k:')
plot(x2_array,fitfun(ci_N(2,1),N_fit.b,x2_array),'k:')
plot (x2_array,x2_array,'r--')
hold off
ylabel('log(N_{detected})')
xlabel('log(N_{catalog})')
axis equal

figure(4)
plot(REW1548_det_all,REW1548_cat_all,'.r'); 
hold on
plot(0:2,0:2,'b--')
xlabel('REW_{detected}'); ylabel('REW_{catalog}')
hold off

figure(5)
plot(REW1548_det_all./REW1550_det_all,REW1548_cat_all./REW1550_cat_all,'.r'); 
hold on
plot(1:2,1:2,'b--')
xlabel('DR_{detected}'); ylabel('DR_{catalog}')
hold off