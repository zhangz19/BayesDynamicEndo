function [x] = updateEta2(x) %, err, updateEtaV
% update delta and gamma separately
global C EV N T1 J  %T

for update_eta = 1:EV.update_etas
    if EV.common_delta == 0 % varying delta
        for i = 1:N
            Lo = zeros(1); Mu = zeros(1,1);
            for t = 2:T1
                Xtild = x.S(t-1,i);
                Lo = Lo + (Xtild'*Xtild)*x.Pre_s;
                Mu = Mu + Xtild'*( x.S(t,i) - C(:,t-1,i)'*x.gamma(:,i) - x.Ome_s*(C(:,t-1,i)-x.theta(:,t,i)) )*x.Pre_s;
            end
            Lo = Lo + diag(1./[x.delta_sigma2]); Mu = Mu + [x.delta_mean]./[x.delta_sigma2];
            Lo = chol(Lo, 'lower');
            Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
            eta = Lo'\Mu;
            x.delta(:,i) = eta(1);
        end
        % for flat prior on delta_mean(varyind):
        %         x.delta_mean = mean(x.delta, 2) + randn(1).*sqrt(x.delta_sigma2/N);
        % % for normal prior with mean 0 and pre_beta
        Mu = sum(x.delta)/x.delta_sigma2;
        Lo = 1/( N/x.delta_sigma2 + EV.pre_eta );
        Mu = Lo*Mu;
        x.delta_mean = Mu + sqrt(Lo)*randn(1);
        bzeros = ( EV.invbeta_delta + 0.5*sum((x.delta - repmat(x.delta_mean, [1,N])).^2, 2) )^-1;
        azeros = EV.alpha_delta + 0.5*N;
        x.delta_sigma2 = 1/gamrnd(azeros, bzeros);
        
    else %common delta
        Lo = zeros(1); Mu = zeros(1,1);
        for i = 1:N
            for t = 2:T1
                Xtild = x.S(t-1,i);
                Lo = Lo + (Xtild'*Xtild)*x.Pre_s;
                Mu = Mu + Xtild'*( x.S(t,i) - C(:,t-1,i)'*x.gamma(:,i) - x.Ome_s*(C(:,t-1,i)-x.theta(:,t,i)) )*x.Pre_s;
            end
        end
        Lo = Lo +  diag(EV.pre_eta).*eye(size(Lo,1));
        Lo = chol(Lo, 'lower');
        Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
        eta = Lo'\Mu;
        x.delta_mean = eta;
        x.delta = repmat(x.delta_mean, [1,N]);
    end
    
    
    if EV.common_gamma == 0 % varying gamma
        for i = 1:N
            Lo = zeros(J); Mu = zeros(J,1);
            for t = 2:T1
                Xtild = C(:,t-1,i)';
                Lo = Lo + (Xtild'*Xtild)*x.Pre_s;
                Mu = Mu + Xtild'*( x.S(t,i) - x.S(t-1,i)*x.delta(i) - x.Ome_s*(C(:,t-1,i)-x.theta(:,t,i)) )*x.Pre_s;
            end
            Lo = Lo + diag(1./[x.gamma_sigma2]); Mu = Mu + [x.gamma_mean]./[x.gamma_sigma2];
            Lo = chol(Lo, 'lower');
            Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
            eta = Lo'\Mu;
            x.gamma(:,i) = eta;
        end
        % for flat prior
        %         x.gamma_mean = mean(x.gamma, 2) + randn([J,1]).*sqrt(x.gamma_sigma2/N);
        % for informative, conjugate prior
        Mu = sum(x.gamma, 2)./x.gamma_sigma2;
        Lo = diag(N./x.gamma_sigma2) + EV.pre_eta*eye(J);
        Lo = chol(Lo, 'lower');
        Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
        x.gamma_mean = Lo'\Mu;
        bzeros = ( EV.invbeta_gamma + 0.5*sum((x.gamma - repmat(x.gamma_mean, [1,N])).^2, 2) ).^-1;
        azeros = EV.alpha_gamma + 0.5*N;
        for j = 1:J
            x.gamma_sigma2(j) = 1./gamrnd(azeros, bzeros(j));
        end
        
    else %common delta
        Lo = zeros(J); Mu = zeros(J,1);
        for i = 1:N
            for t = 2:T1
                Xtild = C(:,t-1,i)';
                Lo = Lo + (Xtild'*Xtild)*x.Pre_s;
                Mu = Mu + Xtild'*( x.S(t,i) - x.S(t-1,i)*x.delta(i) - x.Ome_s*(C(:,t-1,i)-x.theta(:,t,i)) )*x.Pre_s;
            end
        end
        Lo = Lo +  diag(EV.pre_eta).*eye(size(Lo,1));
        Lo = chol(Lo, 'lower');
        Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
        eta = Lo'\Mu;
        x.gamma_mean = eta;
        x.gamma = repmat(x.gamma_mean, [1,N]);
    end
    
    
end

% % step 4.1: or update delta only, update_etas must turn 0
% % or update delta and gamma separately
% % works only for AR model with no gamma terms
% for myendogeneity = 1:(1-EV.endogeneity)
%     %             % note: need to update S for t=0, though no need for t=1~T
%     %             Lo = 1/(delta_mean^2/sigma2_e + pre_S1); % normal prior for S(t=0) with mean 0 and precision pre_S1.
%     %             Mu = Lo*delta_mean*err(1,:)/sigma2_e;
%     %             S0 = sqrt(Lo)*randn(1) + Mu;
%
%     % or if fix S0 = 0
%     S0 = zeros(1, N);
%     tmp = [S0; err];
%
%     if EV.common_delta == 1
%         S1 = reshape(tmp(2:T1,:), [T*N, 1]);
%         S0 = reshape(tmp(1:T,:), [T*N, 1]);
%         Lo = 1/(S0'*S0/x.sigma2_e + EV.pre_eta);
%         Mu = Lo*(S1'*S0)/x.sigma2_e;
%         x.delta_mean = sqrt(Lo)*randn(1) + Mu;
%
%         if EV.fixdelta0 == 1
%             x.delta_mean = 0;
%         end
%
%         x.delta = x.delta_mean*ones(size(x.delta));
%
%     else % update delta as random effects
%         if updateEtaV == 1
%             for i = 1:N
%                 S1 = tmp(2:T1, i);
%                 S0 = tmp(1:T, i);
%                 Lo = 1/(S0'*S0/x.sigma2_e + EV.pre_eta);
%                 Mu = Lo*(S1'*S0)/x.sigma2_e;
%                 x.delta(i) = sqrt(Lo)*randn(1) + Mu;
%             end
%         end
%         x.delta_mean = mean(x.delta, 2) + randn(1).*sqrt(x.delta_sigma2/N);
%         bzeros = 1/( EV.invbeta_delta + 0.5*sum((x.delta - repmat(x.delta_mean, [1,N])).^2, 2) );
%         for j = 1:1
%             x.delta_sigma2(j) = 1./gamrnd(EV.alpha_delta, bzeros);
%         end
%
%     end
% end

end
