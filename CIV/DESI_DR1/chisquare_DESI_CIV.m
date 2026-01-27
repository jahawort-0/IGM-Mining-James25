difference = (this_mu - this_flux).^2
weighteddiff = difference./ (this_mu.^2)
reducedChi2 = mean(weighteddiff);

fit_chi2(this_quasar_ind) = reducedChi2;