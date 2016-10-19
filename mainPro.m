function [] = mainPro(ID)
% A Bayesian dynamic model with semi-varying coefficients
% Zhen Zhang, PhD, zhangquake@gmail.com
% mcc -R -nodisplay -R -singleCompThread -m mainPro.m -v
global Y X C pv pc comind varyind N T T1 J J1 indS indT EV

outofsampletest = 0;
verbose = 0;
nChain = 5; %8 if outofsampletest==1
nVar = 10; %nvar indicates cases here
simu = 0;
simulatedata = 0;

tot = 2e4; burnin = 18e3;

EV.usesubset = 0;
EV.computeDIC = 0;
EV.common_beta = 0; %obsolete
EV.common_eta = 0;
EV.common_delta = 0; % effective only when updateEta2
EV.common_gamma = 0; % effective only when updateEta2
EV.common_phi = 0;

EV.endogeneity = 1;
EV.flatprior4IG = 0;
EV.update_Ses = 1;
EV.update_thetas = 1;
EV.update_etas = 1; % update gamma and delta?
EV.update_betas = 1;
EV.update_sigma2s = 1;
EV.update_phis = 1;
EV.updateMats = 1;
EV.pre_beta = 1e-6;
EV.pre_eta = 1e-6; %EV.pre_delta = 1e-6;  EV.pre_gamma = 1e-6;
EV.pre_phi = 1e-6;
EV.pre_S1 = 1e-6; EV.pre_theta1 = 1e-6;
EV.startT = 1;
EV.singlesite = 0;
% EV.fixdelta0 = 1;

ch = str2double(num2str(ID));
ch0 = ch;
ncase = ceil(ch/(nChain*nVar)); %nvar indicates model id
ch = ch - (ncase-1)*(nChain*nVar);
nvar = ceil(ch/nChain);
ch = ch - (nvar-1)*nChain;
fprintf('scenario = %d, model = %d, chain = %d:\n', [ncase,nvar,ch])

switch nvar
    case 1;
    case 2; EV.common_delta = 1;
    case 3; EV.common_gamma = 1;
    case 4; EV.common_beta = 1;
    case 5; EV.common_phi = 1;
    case 6; EV.endogeneity = 0;
    case 7; EV.endogeneity = 0; EV.common_delta = 1;
    case 8; EV.endogeneity = 0; EV.common_gamma = 1;
    case 9; EV.endogeneity = 0; EV.common_beta = 1;
    case 10; EV.endogeneity = 0; EV.common_phi = 1;
end

load('datFirm.mat')
if outofsampletest == 1
    T = 20; ind0 = find(indT<=T); ind1 = 1:numel(indT); ind1(ind0) = [];
    C1 = C(:,(T+1):end,:); X1 = X(ind1,:); Y1 = Y(ind1,:); indS1 = indS(ind1); indT1 = indT(ind1);
    Y1 = reshape(Y1, [size(C1,2), size(C1,3)]);
    C = C(:,1:T,:); X = X(ind0,:); Y = Y(ind0,:); indS = indS(ind0); indT = indT(ind0);
end

if simulatedata == 1
    rng('default'); rng(12);
    N = 100; T = 24; indS = kron(1:N, ones(1,T))'; NT = N*T;
    X = [normrnd(0,1,[NT,1]), normrnd(3,5,[NT,1])]; C = unifrnd(-1,1,[3, T,N]); Y = normrnd(0,1,[NT,1]);
end

if EV.singlesite == 1
    usesite = ncase;
    X = X(indS==usesite, :); Y = Y(indS==usesite);
    C = C(:, :, usesite);
    indS = indS(indS==usesite); indS = ones(numel(indS),1);
end

N = max(indS); T = length(Y)/N;

p = size(X,2);

% discard quadratic terms in C
J = 4; %size(C,1);
C = C(1:J, :, :);
if outofsampletest == 1
    C1 = C1(1:J,:,:);
end

Y = reshape(Y, [T,N]);
if EV.usesubset == 1
    N = 1; indS = indS(indS<=N);
    Y = Y(:,1:N); T = size(Y,1); NT = N*T; X = X(1:NT,:); C = C(:,:,1:N);
end
NT = N*T; J1 = J+1; T1 = T+1;

nam = strcat('out_', num2str(ncase),'_', num2str(nvar),'_',num2str(ch),'.mat');

%============================ set hyperparameters for priors: START
mean_IG = 0.01; var_IG = 10^4;  alpha_IG = 2+mean_IG^2/var_IG;
invbeta_IG = mean_IG*(alpha_IG-1);
EV.alpha_e = alpha_IG + 0.5*NT;
EV.invbeta_e = invbeta_IG;
EV.alpha_w = alpha_IG + 0.5*NT;
EV.invbeta_w = invbeta_IG;
EV.alpha_beta = alpha_IG + 0.5*N;
EV.invbeta_beta = invbeta_IG;
EV.alpha_gamma = alpha_IG; EV.invbeta_gamma = invbeta_IG;
EV.alpha_delta = alpha_IG + 0.5*N;
EV.invbeta_delta = invbeta_IG;
EV.alpha_phi = alpha_IG; EV.invbeta_phi = invbeta_IG;
if EV.flatprior4IG == 1 % better not for beta
    EV.alpha_e = 0.5*NT; EV.invbeta_e = 0;
    EV.alpha_w = 0.5*NT; EV.invbeta_w = 0;
end
EV.nu_Omega = J1+2   +  NT; % add posterior part
EV.A_Omega = eye(J1);
EV.nu_Psi = J+2  +  NT ; % add posterior part
EV.A_Psi = eye(J);
%============================= set hyperparameters for priors: END




% initialize
% ======================================================================================
if simu == 1
    rng('default'); rng(25);
end
x.beta_mean = 5*ones(p,1);
% beta_sigma2 = .01*ones(p, 1);
x.gamma_mean = 0.5*ones(J,1);
x.gamma_sigma2 = .01*ones(J, 1)*(EV.common_gamma==0);
x.sigma2_e = 1;
if simu == 0
    if EV.singlesite == 0
        x.beta_mean = init_beta;
        x.gamma_mean = init_gamma;
        x.sigma2_e = init_sigmae^2;
        x.beta = init_beta_i;
        x.gamma = init_gamma_i;
        x.beta_sigma2 = init_beta_sigma2;
        x.gamma_sigma2 = init_gamma_sigma2;
        x.delta_mean = init_delta;
        x.beta_sigma2 = init_beta_sigma2;
    else
        x.beta_mean = init_beta_i(:, ncase);
        x.gamma_mean = 0*x.gamma_mean;
        x.sigma2_e = init_sigmae_i(ncase);
    end
    
    if EV.update_etas == 1 % this means we would separate beta and gamma
        varyind = find(varyind == 1);
        pv = numel(varyind); comind = 1:p; comind(varyind) = []; pc = numel(comind);
        
    else
        varyind = find(varyind == 1);
        pv = numel(varyind); comind = 1:p; comind(varyind) = []; pc = numel(comind);
    end
    
    
    if EV.common_beta == 1
        varyind = []; % varyind(1:3); %[];
        pv = numel(varyind); comind = 1:p; pc = p;
    end
end

% set random seeds
rng('default'); rng(ch0*44);%210
x.delta_sigma2 = .01*ones(1, 1)*(EV.common_eta==0);

x.delta = repmat(x.delta_mean, [1,N]);
if EV.common_delta ==0
    x.delta = x.delta + normrnd(0,1,size(x.delta)).*repmat(sqrt(x.delta_sigma2), [1,N]);
end
x.phi_mean = .9*ones(J,1); x.phi_sigma2 = .01*ones(J, 1)*(EV.common_phi==0);
x.phi = repmat(x.phi_mean, [1,N]); x.phi = x.phi + normrnd(0,1,size(x.phi)).*repmat(sqrt(x.phi_sigma2), [1,N]);
x.theta = 10*ones(J,T1,N);
if simu ~= 1
    for t = 2:T1
        for i = 1:N
            for j = 1:J
                x.theta(j,t,i) = C(j,t-1,i)*.9;
            end
        end
    end
    x.theta(:,1,:) = mean(C, 2);
    % sum(sum(sum(theta(:,2:end,:)-C)))
end
x.S = zeros(T1,N);
x.Omega = 1*(eye(J1) + .2); x.Psi = 1*(eye(J) + .2); %eye(J); +.2

Pre = x.Omega\eye(J1);
x.Pre_s = Pre(1,1); % this should = inv(Omega(1,1)-Omega(1,2:J1)*inv(Omega(2:J1,2:J1))*Omega(2:J1,1))
x.Pre_c = Pre(2:J1, 2:J1); % this should = inv(Omega(2:J1,2:J1)-Omega(2:J1,1)*inv(Omega(1,1))*Omega(1,2:J1))
x.Ome_s = -Pre(1,2:J1)/Pre(1,1); %should = Omega(1,2:J1)*inv(Omega(2:J1,2:J1))
x.Ome_c = x.Omega(2:J1,1)/x.Omega(1,1);
x.invPsi = x.Psi\eye(J);

initvec = [...
    reshape(x.beta, [1, numel(x.beta)]), reshape(x.beta_mean, [1, numel(x.beta_mean)]), x.beta_sigma2', x.sigma2_e,...
    reshape(x.delta, [1, numel(x.delta)]), reshape(x.delta_mean, [1, numel(x.delta_mean)]), x.delta_sigma2', ...
    reshape(x.gamma, [1, numel(x.gamma)]), reshape(x.gamma_mean, [1, numel(x.gamma_mean)]), x.gamma_sigma2', ...
    reshape(x.phi, [1, numel(x.phi)]), reshape(x.phi_mean, [1, numel(x.phi_mean)]), x.phi_sigma2',...
    x.Omega(~~tril(x.Omega+5))', x.Psi(~~tril(x.Psi+5))'...
    ];


x = getInits(x, EV);

if simu == 1 % for simulation, use the parameter values above, simulate Y's
    w = zeros(J1, T, N); ksi = zeros(J, T, N); eps = normrnd(0,sqrt(x.sigma2_e),[T,N]);
    err = zeros(T,N);
    for i = 1:N
        err(:,i) = X(indS==i,:)*x.beta(:, i); % err: response surface
        for t = 1:T
            w(:,t,i) = chol(x.Omega,'lower')*normrnd(0,1,[J1,1]);
            ksi(:,t,i) = chol(x.Psi,'lower')*normrnd(0,1,[J,1]);
            x.theta(:,t+1,i) = x.phi(:,i).*x.theta(:,t,i) + ksi(:,t,i);
            C(:,t,i) = x.theta(:,t+1,i) + w(2:J1,t,i);
            x.S(t+1,i) = x.delta(i)*x.S(t,i) + C(:,t,i)'*x.gamma(:,i) + w(1,t,i);
            Y(t,i) = err(t,i) + x.S(t+1,i) + eps(t,i);
        end
        err(:,i) = Y(:,i) - err(:,i); % err: residual Y - Xbeta
    end
else
    % w = normrnd(0, sqrt(Omega(1,1)), [T, N]);
    w = normrnd(0, sqrt(x.sigma2_e), [T, N]);
    if EV.update_Ses == 1
        for i = 1:N
            x.S(1,:) = x.S(1,:) + .01;
            for t = 1:T
                x.S(t+1,i) = x.delta(i)*x.S(t,i)  + w(t,i)+ C(:,t,i)'*x.gamma(:,i);
            end
        end
    end
end
inittheta = reshape(x.theta, [1,numel(x.theta)]);
initS = reshape(x.S, [1,numel(x.S)]);
if EV.update_betas == 0
    err = Y - x.S(2:end,:);
    for i = 1:N
        err(:,i) = err(:,i) - X(indS==i,comind)*x.beta_mean(comind) - X(indS==i,varyind)*x.beta(varyind, i);
    end
end
% ======================================================================================


% MCMC running: store the results, pre-allocation
npara = numel(initvec); ntheta = numel(inittheta); nS = numel(initS);
matPara = nan((tot-burnin), npara);
matTheta = nan((tot-burnin), ntheta);
matS = nan((tot-burnin), nS);
L0s = zeros(tot-burnin, 3); %likelihood, MSE, RMAD for fitted data
Es = [];
if EV.computeDIC == 1
    Es = zeros(tot-burnin, 2);
end
Yhat = [];
if outofsampletest == 1
    Yhat = nan((tot-burnin), N, size(C1,2));
else
    Yhat = nan((tot-burnin), N, size(C,2));
end

checkpoint = 1;
iter0 = 0; t0 = 0; completed = 0;

tic
for iter = (iter0+1):tot
    if verbose == 1
        fprintf('%6d', iter)
    end
    
    %     %step 1: update beta components
    [x, err] = updateBeta(x, 1);
    
    % step 2: update the unexplained variation
    [x] = updateSigma(x, err);
    
    % step 5: update Omega and Psi
    [x] = updateMatrix(x);
    
    % step 3: update S: T by N matrix
    for update_S = 1:EV.update_Ses
        err = err + x.S(2:end,:); % extract old S effect. Here err must be output from updateBeta
        for i = 1:N
            for t = EV.startT:T1
                % t = T1+1-t_foo;
                Lo = 0; Mu = 0;
                if t >1
                    Delta2 = x.delta(i)*x.S(t-1,i) + C(:,t-1,i)'*x.gamma(:,i) + x.Ome_s*(C(:,t-1,i)-x.theta(:,t,i));
                    Lo = Lo + (x.Pre_s + 1/x.sigma2_e); Mu = Mu + (err(t-1,i)/x.sigma2_e + Delta2*x.Pre_s);
                end
                if t < T1
                    Delta3 = x.S(t+1,i) - C(:,t,i)'*x.gamma(:,i) - x.Ome_s*(C(:,t,i)-x.theta(:,t+1,i));
                    Lo = Lo + x.delta(i)^2*x.Pre_s; Mu = Mu + x.delta(i)*x.Pre_s*Delta3;
                    if t == 1
                        Lo = Lo + EV.pre_S1; % normal prior for S(t=0) with mean 0 and precision pre_S1.
                    end
                end
                Lo = 1/Lo;
                x.S(t,i) = Lo*Mu + randn(1)*sqrt(Lo);
            end
        end
        err = err - x.S(2:end,:); % add the new S effect
    end
    
    % step 4: updat delta and gamma components
    [x] = updateEta2(x);  %, err, 1
    
    %     % after updating delta, update S here for prediction later: for AR(1) only
    %     tmp = sqrt(x.sigma2_e) * randn([T, N]);
    %     for t = 2:T1
    %         x.S(t,:) = tmp(t-1,:) + x.delta.*err(t-1,:); % note this is sequentially updating %x.S(t-1,:)
    %     end
    
    % step 6: update phi components
    for update_phi = 1:EV.update_phis
        if EV.common_phi == 0
            for i = 1:N
                Lo = zeros(J); Mu = zeros(J,1);
                for t = 2:T1
                    Lo_t = repmat(x.theta(:,t-1,i), [1,J]).*x.invPsi; Mu = Mu + Lo_t*x.theta(:,t,i);
                    Lo = Lo + Lo_t.*repmat(x.theta(:,t-1,i)', [J,1]);
                end
                Lo = Lo + diag(1./x.phi_sigma2); Mu = Mu + x.phi_mean./x.phi_sigma2;
                Lo = chol(Lo, 'lower');
                Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
                x.phi(:, i) = Lo'\Mu;
            end
            x.phi_mean = mean(x.phi, 2) + randn(size(x.phi_mean)).*sqrt(x.phi_sigma2/N);
            bzeros = ( EV.invbeta_phi + 0.5*sum((x.phi - repmat(x.phi_mean, [1,N])).^2, 2) ).^-1;
            azeros = EV.alpha_phi + 0.5*N;
            for j = 1:J
                x.phi_sigma2(j) = 1./gamrnd(azeros, bzeros(j));
            end
        else
            Lo = zeros(J); Mu = zeros(J,1);
            for i = 1:N
                for t = 2:T1
                    Lo_t = repmat(x.theta(:,t-1,i), [1,J]).*x.invPsi; Mu = Mu + Lo_t*x.theta(:,t,i);
                    Lo = Lo + Lo_t.*repmat(x.theta(:,t-1,i)', [J,1]);
                end
            end
            Lo = Lo + diag(EV.pre_phi);
            Lo = chol(Lo, 'lower');
            Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
            x.phi_mean = Lo'\Mu;
            x.phi = repmat(x.phi_mean, [1,N]);
        end
    end
    
    % step 7: update theta components
    %tic
    for update_theta = 1:EV.update_thetas
        for i = 1:N
            for t = EV.startT:T1
                Lo = zeros(J); Mu = zeros(J,1);
                if t >1
                    Delta2 = x.S(t,i) - x.delta(i)*x.S(t-1,i) - C(:,t-1,i)'*x.gamma(:,i);
                    Lo = Lo + (x.Pre_c + x.invPsi); Mu = Mu + x.Pre_c*(C(:,t-1,i) - x.Ome_c*Delta2) + x.invPsi*(x.phi(:,i).*x.theta(:,t-1,i));
                end
                if t < T1
                    Lo_t = repmat(x.phi(:,i), [1,J]).*x.invPsi; Mu = Mu + Lo_t*x.theta(:,t+1,i);
                    Lo = Lo + Lo_t.*repmat(x.phi(:,i)', [J,1]);
                    if t == 1 % normal prior for gamma (t=0)
                        Lo = Lo + eye(J)*EV.pre_theta1;
                    end
                end
                Lo = chol(Lo, 'lower'); Mu = Lo\Mu; Mu = Mu + randn(size(Mu));
                x.theta(:,t,i) = Lo'\Mu;
            end
        end
    end
    %toc
    
    % store results
    if iter > burnin
        matPara((iter-burnin),:) = [...
            reshape(x.beta, [1, numel(x.beta)]), reshape(x.beta_mean, [1, numel(x.beta_mean)]), x.beta_sigma2', x.sigma2_e, ...
            reshape(x.delta, [1, numel(x.delta)]), reshape(x.delta_mean, [1, numel(x.delta_mean)]), x.delta_sigma2', ...
            reshape(x.gamma, [1, numel(x.gamma)]), reshape(x.gamma_mean, [1, numel(x.gamma_mean)]), x.gamma_sigma2', ...
            reshape(x.phi, [1, numel(x.phi)]), reshape(x.phi_mean, [1, numel(x.phi_mean)]), x.phi_sigma2',...
            x.Omega(~~tril(x.Omega+5))', x.Psi(~~tril(x.Psi+5))'...
            ];
        
        err2 = reshape(err, [1,numel(err)]);
        L0s(iter-burnin, :) = [ sum( log(normpdf(err2, 0, sqrt(x.sigma2_e))) ), mean(err2.^2), ...
            mean(reshape(abs(err./Y), [1,numel(err)])) ];
        
        if outofsampletest == 1
            Lo = chol(x.invPsi, 'lower');
            % The predictive distribution simulation needs a further MC
            % integration evaluation over the latent process.
            Len = 50;
            for i = 1:N
                Yhat(iter-burnin, i, :) = X1(indS1==i,varyind)*x.beta(varyind, i) + X1(indS1==i,comind)*x.beta_mean(comind) ...
                    + randn([size(Yhat,3),1])*sqrt(x.sigma2_e);
                for subrep = 1:Len
                    theta0 = x.phi(:,i).*x.theta(:,T1,i);
                    S0 = x.delta(i)*x.S(T1,i) + C(:,T,i)'*x.gamma(:,i);
                    for t = 1:size(C1,2)
                        theta = theta0 + Lo'\randn([J,1]);
                        S = S0 + x.Ome_s*(C1(:,t,i)-theta) + randn(1)/sqrt(x.Pre_s);
                        Yhat(iter-burnin, i, t) = Yhat(iter-burnin, i, t) + S/Len;
                        theta0 = x.phi(:,i).*theta;
                        S0 = x.delta(i)*S + C1(:,t,i)'*x.gamma(:,i);
                    end
                end
            end
            tmp = (squeeze(Yhat(iter-burnin, :, :)) - Y1').^2;
            tmp = sqrt([   mean(reshape(err.^2, [1,numel(err)])),  mean(reshape(tmp,[1,numel(tmp)]))   ]); %RMSE
            if verbose==1
                fprintf('Lik = %3.3f, RMSE = %3.3f, out-of-sample RMSE = %3.3f\n', [L0s(iter-burnin, 1),tmp])
            end
            
        else % predict the full Y
            Yhat(iter-burnin, :, :) = (Y - err)';
        end
        
        % for computing DIC
        for mycomputeDIC = 1:EV.computeDIC
            Es(iter-burnin, 1) = L0s(iter-burnin,1);% complete likelihood
            if pv >0
                tmp = ( x.beta(varyind, :) - repmat(x.beta_mean(varyind), [1, N]) )./ repmat(sqrt(x.beta_sigma2(varyind)), [1, N]);
                Es(iter-burnin, 1) = Es(iter-burnin, 1) + sum(sum( log(normpdf(tmp)) )); % add likelihood for beta
                if EV.common_eta == 0
                    tmp = ( x.delta - repmat(x.delta_mean, [1, N]) )./ repmat(sqrt(x.delta_sigma2), [1, N]);
                    Es(iter-burnin, 1) = Es(iter-burnin, 1) + sum(sum( log(normpdf(tmp)) )); % add likelihood for delta
                end
                
                len = 30;
                beta_b = zeros(size(x.beta)); beta_mean_b = zeros(size(x.beta_mean)); beta_sigma2_b = zeros(size(x.beta_sigma2));
                sigma2_e_b = zeros(size(x.sigma2_e));
                delta_b = zeros(size(x.delta)); delta_mean_b = zeros(size(x.delta_mean)); delta_sigma2_b = zeros(size(x.delta_sigma2));
                
                x0 = x;
                for i = 1:len
                    [x, err] = updateBeta(x, 0); beta_b = beta_b+x.beta; beta_mean_b = beta_mean_b+x.beta_mean; beta_sigma2_b = beta_sigma2_b+x.beta_sigma2;
                    [x] = updateSigma(x, err); sigma2_e_b = sigma2_e_b+x.sigma2_e;
                    if EV.common_eta == 1
                        [x] = updateEta(x, err, 1); delta_b = delta_b+x.delta; delta_mean_b = delta_mean_b+x.delta_mean; delta_sigma2_b = delta_sigma2_b+x.delta_sigma2;
                    else % do not update x.delta
                        [x] = updateEta(x, err, 0); delta_b = delta_b+x.delta; delta_mean_b = delta_mean_b+x.delta_mean; delta_sigma2_b = delta_sigma2_b+x.delta_sigma2;
                    end
                end
                x = x0;
                
                beta_b=beta_b/len; beta_mean_b=beta_mean_b/len; beta_sigma2_b=beta_sigma2_b/len;
                delta_b=delta_b/len; delta_mean_b=delta_mean_b/len; delta_sigma2_b=delta_sigma2_b/len;
                sigma2_e_b=sigma2_e_b/len;
                err = zeros(T,N);
                for i = 1:N
                    err(:,i) = Y(:,i) - X(indS==i,varyind)*beta_b(varyind, i) - X(indS==i,comind)*beta_mean_b(comind); % - Smean(2:end,i);
                end
                tmp = err; %(2:end,:);
                % for ordinary AR(1) model only
                tmp = [tmp(1,:); tmp(2:T,:) - repmat(delta_b,[T-1,1]).*tmp(1:(T-1),:)];
                tmp = tmp/sqrt(sigma2_e_b);
                Es(iter-burnin, 2) = sum(sum( log(normpdf(tmp)) ));
                
                tmp = (beta_b(varyind, :) - repmat(beta_mean_b(varyind), [1, N]) )./ repmat(sqrt(beta_sigma2_b(varyind)), [1, N]);
                Es(iter-burnin, 2) = Es(iter-burnin, 2) + sum(sum( log(normpdf(tmp)) )); % add likelihood for beta
                if EV.common_eta == 0
                    tmp = ( delta_b - repmat(delta_mean_b, [1, N]) )./ repmat(sqrt(delta_sigma2_b), [1, N]);
                    Es(iter-burnin, 2) = Es(iter-burnin, 2) + sum(sum( log(normpdf(tmp)) )); % add likelihood for delta
                end
                
            end
        end
        
    end
    
    if toc > 14000 && checkpoint == 1
        iter0 = iter; CPUtime = toc; t0 = t0 + CPUtime/60; xN = x;
        save(nam,'matPara','t0','iter0','completed','initvec','inittheta','initS','L0s','Es','Yhat','EV','xN') %'matTheta','matS',
        checkpoint = 0;
    end
    
end


CPUtime = toc; CPUtime = t0 + CPUtime/60;
completed = 1; xN = x;
fprintf('\n%d iterations are done with elapsed time %.2f minutes.\n', tot, CPUtime)
save(nam,'matPara','t0','iter0','completed','initvec','inittheta','initS','L0s','Es','Yhat','EV','xN') %,'matTheta','matS'
end





