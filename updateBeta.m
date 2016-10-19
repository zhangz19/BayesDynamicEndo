function [x, err] = updateBeta(x, updateBetaV)
% if updateBetaV ==1: not updating beta(varyind,:), but update beta_mean
global EV Y X pv pc comind varyind N indS %T XtX
for update_beta = 1:EV.update_betas
    
    err = Y - x.S(2:end,:); % T by N
    
    % updating varying-beta components
    if pv > 0
        
        if updateBetaV == 1
            for i = 1:N
                Xt = X(indS==i,varyind); %T by p
                %                 % XtX for ordinary AR(1) model
                %                 Xt = [Xt(1,:); Xt(2:T,:) - x.delta(i)*Xt(1:(T-1),:)];
                XtX = Xt'*Xt;
                Lo = XtX/x.sigma2_e + diag(1./x.beta_sigma2(varyind));
                tmp = err(:,i);
                if pc > 0
                    tmp = tmp - X(indS==i,comind)*x.beta_mean(comind);
                end
                %                 % for ordinary AR(1) model only
                %                 tmp = [tmp(1); tmp(2:T) - x.delta(i)*tmp(1:(T-1))];
                % Mu = X(indS==i,varyind)'*tmp/sigma2_e + beta_mean(varyind)./beta_sigma2(varyind);
                Mu = Xt'*tmp/x.sigma2_e + x.beta_mean(varyind)./x.beta_sigma2(varyind);
                Lo = chol(Lo, 'lower');
                Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
                x.beta(varyind, i) = Lo'\Mu;
                err(:, i) = err(:, i) - X(indS==i,varyind)*x.beta(varyind, i);
            end
        end
        
        % % for flat prior on beta_mean(varyind):
        %         x.beta_mean(varyind) = mean(x.beta(varyind,:), 2) + randn([pv,1]).*sqrt(x.beta_sigma2(varyind)/N);
        % % for normal prior with mean 0 and pre_beta
        Mu = sum(x.beta(varyind,:), 2)./x.beta_sigma2(varyind);
        Lo = diag(N./x.beta_sigma2(varyind)) + EV.pre_beta*eye(pv);
        Lo = chol(Lo, 'lower');
        Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
        x.beta_mean(varyind) = Lo'\Mu;
        
        for j1 = 1:pv
            j = varyind(j1);
            x.beta_sigma2(j) = 1./gamrnd(EV.alpha_beta,  1/( EV.invbeta_beta + 0.5*sum((x.beta(j,:) - repmat(x.beta_mean(j), [1,N])).^2, 2) )  );
        end
    end
    
    % updating common-beta components
    if pc > 0
        Lo = zeros(pc); Mu = zeros(pc,1);
        for i = 1:N
            Xt = X(indS==i,comind); %T by p
            %             Xt = [Xt(1,:); Xt(2:T,:) - x.delta(i)*Xt(1:(T-1),:)];
            XtX = Xt'*Xt;
            Lo = Lo + XtX/x.sigma2_e;
            tmp = err(:,i);
            %             % for ordinary AR(1) model only
            %             tmp = [tmp(1); tmp(2:T) - x.delta(i)*tmp(1:(T-1))];
            Mu = Mu + Xt'*tmp/x.sigma2_e;
        end
        Lo = Lo +  EV.pre_beta*eye(pc);
        Lo = chol(Lo, 'lower');
        Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
        x.beta_mean(comind) = Lo'\Mu;
        for i = 1:N
            err(:, i) = err(:, i) - X(indS==i,comind)*x.beta_mean(comind);
        end
        x.beta(comind,:) = repmat(x.beta_mean(comind), [1,N]);
    end
    
end
end
