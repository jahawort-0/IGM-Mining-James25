%If plots not showing up paste into command: set(0, 'DefaultFigureRenderer', 'painters')
for i=1:30 %numel(z_qsos)
    figure(i)

    wave = all_wavelengths{i};
    flux = all_flux{i};
    ivar = 1./all_noise_variance{i};
    tid = all_tid_dr1{i};
    z = z_qsos(i);
    
    wave = wave./(1+z);
    good = flux.*sqrt(ivar)>0;
    wave = wave(good);
    flux = flux(good);

    plot(wave,flux,'k-'); hold on
    xline(1550,'r-')
    xline(1548,'r-')
    xlim([1310,1555])
    xlabel('Wavelength')
    ylabel('Flux')
    title(sprintf('%s, i=%d',tid,i))
end

%strong absorbers in i=7,13,20