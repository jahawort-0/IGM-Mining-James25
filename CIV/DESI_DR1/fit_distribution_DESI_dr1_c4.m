fit_difference = this_flux - this_mu; %difference between flux and continuum
avgflux = mean(this_flux);  %normalization constant
mean_diff = mean(fit_difference) ./ avgflux;
std_diff = std(fit_difference) ./ avgflux;

if (plotting == 1)
    fig = figure('visible', 'off');
    histogram(fit_difference,20)
    xlabel('flux difference from continuum model')
    title(sprintf('AVG = %.3f    std = %.3f',mean_diff, std_diff))
    xline(0,'m-','LineWidth',3)


    fid = sprintf('%s/%s/plt/hist-ind-%d-%s.png', base_directory, releaseTest, all_quasar_ind, all_QSO_ID_dr1{all_quasar_ind});
    exportgraphics(fig, fid,'Resolution', 400,'BackgroundColor','w');
    %clf();
end

all_fit_mean_diff(all_quasar_ind) = mean_diff;
all_fit_std_diff(all_quasar_ind) = std_diff;
