function [B,chi,absorber_gate_mask] = absorber_gate(flux,wave,continuum,doublet,z_abs,sig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% arguments (Input)
%     flux %
%     wave
%     continuum
%     doublet,
%     z_abs
% end
% 
% arguments (Output)
%     B
% end

% define value commands for troubleshooting
% flux = this_flux; wave = this_wavelengths; continuum = this_mu; doublet = c4_muL2;
% z_abs = map_z_c4L2(this_quasar_ind, num_c4); sig = map_sigma_c4L2(this_quasar_ind, num_c4);

%define constants and useful conversion formulas
civ_1548_wavelength = 1548.1949462890625; %constants
civ_1550_wavelength =  1550.77001953125;
speed_of_light = 299792458;
kms_to_z = @(kms) (kms * 1000) / speed_of_light; %functions to convert quantities
wave2z = @(wavelength) (wavelength / civ_1548_wavelength) - 1;
z2wave = @(z) civ_1548_wavelength(z+1);


%define mask by sigma
%sigma_z = sig / speed_of_light; %redshift width of single absorption line
    %define interval we are looking at
    %min is 1548 - 3sigma
    %max is 1550 + 3sigma
%wave_min = civ_1548_wavelength*(1+z_abs-(0.3*sigma_z));
%wave_max = civ_1550_wavelength*(1+z_abs+(0.3*sigma_z));


%OR Define mask by absorption profile difference
absorption = continuum - doublet; %all zero except for two peaks
absorption = absorption./max(absorption); %normalize to max height of 1
d_abs = 0.001; % take all those with abs above this value
abs_range = absorption>=d_abs; %filter abs, this is almost our mask, but we need to find the edges
abs_mask = find(abs_range==1); %find indices of all above d_abs
abs_lower = abs_mask(1)-1; abs_upper = abs_mask(end)+1; %take lowest index-2, highest index+2
wave_min = wave(abs_lower); wave_max = wave(abs_upper); %find wave array value

mask = (wave>wave_min) & (wave<wave_max); absorber_gate_mask = mask; %make mask
wave_mask = wave(mask); %apply mask
flux_mask = flux(mask);
continuum_mask = continuum(mask);
doublet_mask = doublet(mask);

%plotting for troubleshooting
% fig = figure('visible', 'on');
% plot(wave_mask,flux_mask,'b-'); hold on
% plot(wave_mask,continuum_mask,'m-')
% plot(wave_mask,doublet_mask,'r-'); hold off
% plot(wave_mask,flux_mask-doublet_mask); hold on
% plot(wave_mask,doublet_sub); hold off

doublet_sub = continuum_mask - doublet_mask; %area between continuum and doublet
flux_sub = continuum_mask - flux_mask;
B = sum(abs(flux_mask-doublet_mask)) /sum(doublet_sub); %divide difference between data and model by expected area
% B = 0   the flux exactly matches the doublet
% B = 1   The area between flux and doublet is equal to area between
% continuum and doublet (ex. flux~continuum)
% B > 1   The flux is far removed from the doublet


chi = sum((flux_mask-doublet_mask).^2)/sum((doublet_sub.^2));
%chi squared value is similar to B, but weights larger deviations more
%harshly (x^2 > |x|)

end