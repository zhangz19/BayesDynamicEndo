function [x] = updateSigma(x, err)
global EV %T
for update_sigma2 = 1:EV.update_sigma2s
    tmp = err;
    %     % for ordinary AR(1) model only
    %     tmp = [tmp(1,:); tmp(2:T,:) - repmat(x.delta,[T-1,1]).*tmp(1:(T-1),:)];
    x.sigma2_e = 1./gamrnd( EV.alpha_e,  1/( EV.invbeta_e + 0.5*sum(sum(tmp.^2)) ) );
end
end
