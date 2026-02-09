% Load Data
load('/home/cosmic/IGM-Mining-James/output/dr1/processed/processed_sigma_125_Spline_2200-2300_inject_abs.mat')
load('/home/cosmic/IGM-Mining-James/output/dr1/processed/preloaded_qsos_inject_abs.mat')
load('/home/cosmic/IGM-Mining-James/output/dr1/processed/SALSA_catalog.mat')
% load('/home/cosmic/IGM-Mining-Aaron/MGII/output/dr1/processed/catalog.fits')

%% Z catalog vs. Z detected

valid_abs = masked_abs_pixels == 0;
detected_abs = (p_c4 > 0.85);
detected_abs_z = map_z_c4L2.*detected_abs;
% z = 1.5
valid_abs_15 = valid_abs(:,1);
all_checks_15 = abs(detected_abs_z - 1.5) < 0.1;
detected_abs_z_15 = detected_abs_z.*all_checks_15;
detected_abs_z_15(:,1:max_c4) = detected_abs_z_15(:,1:10).*valid_abs_15;
[row_15,col_15] = find(detected_abs_z_15>0);
ind_15 = find(detected_abs_z_15>0);
% z = 2.0
valid_abs_20 = valid_abs(:,2);
all_checks_20 = abs(detected_abs_z - 2.0) < 0.1;
detected_abs_z_20 = detected_abs_z.*all_checks_20;
detected_abs_z_20(:,1:max_c4) = detected_abs_z_20(:,1:10).*valid_abs_20;
[row_20,col_20] = find(detected_abs_z_20>0);
ind_20 = find(detected_abs_z_20>0);
% z = 3.0
valid_abs_30 = valid_abs(:,3);
all_checks_30 = abs(detected_abs_z - 3.0) < 0.1;
detected_abs_z_30 = detected_abs_z.*all_checks_30;
detected_abs_z_30(:,1:max_c4) = detected_abs_z_30(:,1:10).*valid_abs_30;
[row_30,col_30] = find(detected_abs_z_30>0);
ind_30 = find(detected_abs_z_30>0);
% z = 4.0
valid_abs_40 = valid_abs(:,4);
all_checks_40 = abs(detected_abs_z - 4.0) < 0.1;
detected_abs_z_40 = detected_abs_z.*all_checks_40;
detected_abs_z_40(:,1:max_c4) = detected_abs_z_40(:,1:10).*valid_abs_40;
[row_40,col_40] = find(detected_abs_z_40>0);
ind_40 = find(detected_abs_z_40>0);
% z = 5.0
valid_abs_50 = valid_abs(:,5);
all_checks_50 = abs(detected_abs_z - 5.0) < 0.1;
detected_abs_z_50 = detected_abs_z.*all_checks_50;
detected_abs_z_50(:,1:max_c4) = detected_abs_z_50(:,1:10).*valid_abs_50;
[row_50,col_50] = find(detected_abs_z_50>0);
ind_50 = find(detected_abs_z_50>0);

z_detected = [detected_abs_z_15(ind_15);detected_abs_z_20(ind_20);...
    detected_abs_z_30(ind_30);detected_abs_z_40(ind_40);...
    detected_abs_z_50(ind_50)];
z_cat = [sim_z_15(row_15);sim_z_20(row_20);sim_z_30(row_30);...
    sim_z_40(row_40);sim_z_50(row_50)];
y = @(a,x)a*x;
z_fit = fit(z_detected,z_cat,y);
figure(1)
plot(z_detected,z_cat,'.r')
hold on
plot(linspace(0,2.2,1000),y(z_fit.a,linspace(0,2.2,1000)),'k--')
xlabel('z_{detected}')
ylabel('z_{catalog}')