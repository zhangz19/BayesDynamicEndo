function [] = sumPro()
% summary function for
% mixed-effects dynamic model with varying slopes
% global ys  W M E indT varyind indS

continueIt = 1;
saveInits = 0;
nvars = 1:10; %[1,2,3,4];
ncases = 1;% [1:30]; %[1,2,3,4]; % vague
% chs =  [1:5];
chs =  [1,5,8]; % for pred run
% chs= [2,6,7];
% usemodel = [1,2,3];
usemodel = [1,5,6]; % for pred run
usesubset = 0;
niter =  2e3;
n0 = niter;
burn = 0;
thin = 1;
datapath = './';
figpath = './';
monitor_paras = 0;
monitor_thetas = 0;
monitor_Ses = 0;
singlesite = 0;
common_beta = 0;
fixdelta0 = 1;
computeDIC = 0;
outofsampletest = 1;
nch = numel(chs);
nsample = (niter-burn)/thin;
tot = nch*nsample;

disp(chs)

for ncas = ncases
    Measures = nan(numel(nvars), 3+2);
    paras = cell(1,numel(usemodel));
    Ys = cell(1,numel(usemodel));
    ysloop = 1;
    
    for nvar = nvars
        
        fprintf('model = %3d\n',nvar);
        filename = strcat(figpath,'fig',num2str(ncas),'_',num2str(nvar),'.ps');
        
        simu = 0;
        simulatedata = 0;
        if ncas == 3
            simu = 0; simulatedata = 0;
        end
        if ncas == 4
            simu = 0; simulatedata = 0;
        end
        
        load('datFirm.mat')
        Y1 = [];
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
        if singlesite == 1
            usesite = ncas;
            X = X(indS==usesite, :); Y = Y(indS==usesite);
            C = C(:, :, usesite);
            indS = indS(indS==usesite); indS = ones(numel(indS),1);
        end
        
        N = max(indS); T = length(Y)/N; p = size(X,2);
        
        % J = size(C,1);
        J = 4; C = C(1:J, :, :);
        
        varyind = find(varyind==1);
        if common_beta == 1
            varyind = [];
        end
        pv = numel(varyind); comind = 1:p; comind(varyind) = []; pc = numel(comind);
        Npv = N*pv;
        
        Y = reshape(Y, [T,N]);
        if usesubset == 1
            N = 1; indS = indS(indS<=N);
            Y = Y(:,1:N); T = size(Y,1); NT = N*T; X = X(1:NT,:); C = C(:,:,1:N);
        end
        NT = N*T; J1 = J+1; T1 = T+1; Np = N*p; NJ = N*J;
        
        % check the convergence
        for namdefine = 1:1
            Omega = eye(J1); nO = numel(Omega(~~tril(Omega+5)));
            Psi = eye(J); nP = numel(Psi(~~tril(Psi+5)));
            
            nams = [  strcat(repmat({'\beta (k='},[1,Np]),cellfun(@num2str,num2cell(kron(ones(1,N),1:p)),'UniformOutput', false), ...
                repmat({', i='},[1,Np]), cellfun(@num2str,num2cell(kron(1:N, ones(1,p))),'UniformOutput', false), repmat({')'},[1,Np])), ...
                strcat(repmat({'\beta_m (k='},[1,p]),cellfun(@num2str,num2cell(kron(ones(1,1),1:p)),'UniformOutput', false), repmat({')'},[1,p])), ...
                strcat(repmat({'\beta_s (k='},[1,p]),cellfun(@num2str,num2cell(kron(ones(1,1),1:p)),'UniformOutput', false), repmat({')'},[1,p])), ...
                '\sigma2_e', ...
                strcat(repmat({'\delta (k='},[1,N]),cellfun(@num2str,num2cell(kron(ones(1,N),1)),'UniformOutput', false), ...
                repmat({', i='},[1,N]), cellfun(@num2str,num2cell(kron(1:N, ones(1,1))),'UniformOutput', false), repmat({')'},[1,N])), ...
                strcat(repmat({'\delta_m (k='},[1,1]),cellfun(@num2str,num2cell(kron(ones(1,1),1:1)),'UniformOutput', false), repmat({')'},[1,1])), ...
                strcat(repmat({'\delta_s (k='},[1,1]),cellfun(@num2str,num2cell(kron(ones(1,1),1:1)),'UniformOutput', false), repmat({')'},[1,1])), ...
                strcat(repmat({'\gamma (k='},[1,NJ]),cellfun(@num2str,num2cell(kron(ones(1,N),1:J)),'UniformOutput', false), ...
                repmat({', i='},[1,NJ]), cellfun(@num2str,num2cell(kron(1:N, ones(1,J))),'UniformOutput', false), repmat({')'},[1,NJ])), ...
                strcat(repmat({'\gamma_m (k='},[1,J]),cellfun(@num2str,num2cell(kron(ones(1,1),1:J)),'UniformOutput', false), repmat({')'},[1,J])), ...
                strcat(repmat({'\gamma_s (k='},[1,J]),cellfun(@num2str,num2cell(kron(ones(1,1),1:J)),'UniformOutput', false), repmat({')'},[1,J])), ...
                strcat(repmat({'\phi (k='},[1,NJ]),cellfun(@num2str,num2cell(kron(ones(1,N),1:J)),'UniformOutput', false), ...
                repmat({', i='},[1,NJ]), cellfun(@num2str,num2cell(kron(1:N, ones(1,J))),'UniformOutput', false), repmat({')'},[1,NJ])), ...
                strcat(repmat({'\phi_m (k='},[1,J]),cellfun(@num2str,num2cell(kron(ones(1,1),1:J)),'UniformOutput', false), repmat({')'},[1,J])), ...
                strcat(repmat({'\phi_s (k='},[1,J]),cellfun(@num2str,num2cell(kron(ones(1,1),1:J)),'UniformOutput', false), repmat({')'},[1,J])), ...
                strcat(repmat({'\Omega_{'},[1,nO]),cellfun(@num2str,num2cell(1:nO),'UniformOutput', false), repmat({'}'},[1,nO])), ...
                strcat(repmat({'\Psi_{'},[1,nP]),cellfun(@num2str,num2cell(1:nP),'UniformOutput', false), repmat({'}'},[1,nP])), ...
                ];
            
            namsind = [  1*ones(1,Np), ...
                2*ones(1,2*p), ...
                3, ... %'\sigma2_e', ...
                4*ones(1,N),...%strcat(repmat({'\delta (k='},[1,N]),cellfun(@num2str,num2cell(kron(ones(1,N),1)),'UniformOutput', false), repmat({', i='},[1,N]), cellfun(@num2str,num2cell(kron(1:N, ones(1,1))),'UniformOutput', false), repmat({')'},[1,N])), ...
                5,...%strcat(repmat({'\delta_m (k='},[1,1]),cellfun(@num2str,num2cell(kron(ones(1,1),1:1)),'UniformOutput', false), repmat({')'},[1,1])), ...
                6,...%strcat(repmat({'\delta_s (k='},[1,1]),cellfun(@num2str,num2cell(kron(ones(1,1),1:1)),'UniformOutput', false), repmat({')'},[1,1])), ...
                7*ones(1,NJ), ... % for gamma
                8*ones(1,2*J), ...
                11*ones(1,NJ), ... % for phi
                12*ones(1,2*J), ...
                9*ones(1,nO), ...%strcat(repmat({'\Omega_{'},[1,nO]),cellfun(@num2str,num2cell(1:nO),'UniformOutput', false), repmat({'}'},[1,nO])), ...
                10*ones(1,nP), ...strcat(repmat({'\Psi_{'},[1,nP]),cellfun(@num2str,num2cell(1:nP),'UniformOutput', false), repmat({'}'},[1,nP])), ...
                ];
            %sum(namsind ~=1 & namsind ~=4 & namsind ~=7 & namsind ~= 11)
            
        end
        npara1 = length(nams);
        npara2 = J*T1*N;
        npara3 = T1*N;
        
        matParas = nan(nch, npara1, nsample);
        matThetas = nan(nch, npara2, nsample);
        matSes = nan(nch, npara3, nsample);
        
        L0sAll = nan(tot, 3);
        EsAll = nan(tot, 2);
        YhatAll = nan(tot, N, size(C1,2));
        
        % combine chains
        usesample = (burn+1):thin:n0;
        for ch = 1:nch
            load(strcat(datapath,'out_', num2str(ncas),'_', num2str(nvar),'_',num2str(chs(ch)),'.mat'))
            tmp = matPara(usesample,:)';
            if sum(sum(sum(isnan(tmp))))>0
                fprintf('scenario = %d, model = %d, chain = %d:\n', [ncas,nvar,ch])
            else
                for j = 1:npara1
                    matParas(ch, j, :) = matPara(usesample,j)';
                end
                %             for j = 1:npara2
                %                 matThetas(ch, j, :) = matTheta(usesample,j)';
                %             end
                %             for j = 1:npara3
                %                 matSes(ch, j, :) = matS(usesample,j)';
                %             end
                L0sAll((ch-1)*nsample + (1:nsample), :) = L0s(usesample, :);
                %             if computeDIC == 1
                %                 EsAll((ch-1)*nsample + (1:nsample),:) = Es(usesample, :);
                %             end
                YhatAll((ch-1)*nsample + (1:nsample), :,:) = Yhat(usesample,:, :);
            end
        end
        %         matParas(isnan(matParas)) = [];
        %         L0sAll(isnan(matParas)) = [];
        %         YhatAll(isnan(matParas)) = [];
        
        YhatAll = permute(YhatAll, [3,2,1]);
        err = repmat(Y1, [1,1,size(YhatAll,3)]) - YhatAll;
        
        %         tmp = zeros(1,nch);
        %         for ch = 1:nch
        %             vec = err(:,:,(ch-1)*nsample+(1:nsample)); vec = reshape(vec, [1, numel(vec)]);
        %             tmp(ch) = mean(vec.^2);
        %         end
        %         sqrt(mean(tmp([1,2,3,5,6])))
        
        err = reshape(err,  [1, numel(err)]);
        err2 = abs( (repmat(Y1, [1,1,size(YhatAll,3)]) - YhatAll)./repmat(Y1, [1,1,size(YhatAll,3)]) );
        %         [~, I] = sort(squeeze( sum(sum(err2.^2, 1),2) )) ;
        %         inds = I; %(1:4e3);
        err2 = reshape(err2,  [1, numel(err2)]);
        Measures(nvar,:) = [mean(-2*L0sAll(:,1)), sqrt(mean(L0sAll(:,2))), mean(L0sAll(:,3)), sqrt(mean(err(:).^2)), mean(err2)];
        
        if any(usemodel == nvar)
            Ys{ysloop} = YhatAll;
        end
        
        if continueIt == 1 %&& nvar == nvars(1)
            % plot(L0sAll')
            EsAll = mean(EsAll, 1);
            
            nbin = 36;
            
            K0 = 1;
            % investigating parameters
            for monitor_para = 1:monitor_paras % turn 1:1 to activate
                vars = sum(var(matParas, [], 3),1); inds = find(vars > 1e-20);
                %npara1 = length(inds); %inds = 1:numel(vars); %
                %matParas = matParas(:, inds,:); initvec = initvec(inds); nams = nams(inds);
                init_i = 1;
                %                 if common_beta == 1
                %                     init_i = 1;
                %                 end
                for j = init_i:size(matParas,2) %944
                    if namsind(j) > 1 && any(inds==j)
                        vec = reshape(matParas(:,j,:),[nch, nsample])';
                        subplot(6,6,K0), plot(vec)
                        title(nams(j),'FontSize',8)
                        %                 if simu == 1
                        %h = refline(0, initvec(j)); set(h, 'LineWidth',2,'Color','k')
                        ylim([min([min(vec),initvec(j)-1e-4]), max([max(vec), initvec(j)+1e-4])])
                        %                 end
                        axis tight; %axis off;
                        K0 = K0 + 1;
                        orient landscape
                        if j == nbin
                            print('-painters', '-dpsc2', '-r600', filename)
                        end
                        if K0 > nbin
                            print('-painters', '-dpsc2', '-r600', '-append', filename)
                            K0 = 1;
                            for ii = 1:nbin % set break point here for monitoring parameters by bins
                                subplot(6,6,ii); plot(nan)
                            end
                        end
                    end
                end
                
            end
            
            R1 = psrf(matParas); %disp(R1)
            Paras = nan(nch*nsample, npara1);
            for j = 1:npara1
                Paras(:,j) = reshape(matParas(:,j,:),[nch*nsample 1]);
            end
            
            %             % for AR(1): need S
            %             Ses = nan(nch*nsample, npara3);
            %             for j = 1:npara3
            %                 Ses(:,j) = reshape(matSes(:,j,:),[nch*nsample 1]);
            %             end
            %             Smean = reshape(mean(Ses, 1), [T1, N]);
            
            mat = zeros(npara1, 3);
            
            for i = 1:npara1
                mat(i,1) = mean(Paras(:,i));
                %     mat(i,2) = std(Paras(:,i));
                %     [lb, ub] = FindHPDset(Paras(:,i), 0.95, []);
                %     if isempty(lb)
                %         disp('here')
                %     end
                %     mat(i,3) = lb; mat(i,4) = ub;
                mat(i,3) = quantile(Paras(:,i),.975);
                mat(i,2) = quantile(Paras(:,i),.025);
                % mat(i,5) = quantile(Paras(:,i),.5);
            end
            %tab1 = [char(nams), num2str([initvec', mat, (initvec'>mat(:,2)).*(initvec'<mat(:,3)) ], 5)];
            %disp(tab1)
            
            if any(usemodel == nvar)
                paras{ysloop} = mat;
                ysloop = ysloop+1;
            end
            
            % save step results as intial values
            if saveInits == 1
                initvec = mat(:,1)';
                save('Inits.mat', 'initvec','namsind')
            end
            
            % investigating latent variables theta
            for monitor_theta = 1:monitor_thetas % turn 1:1 to activate
                R2 = psrf(matThetas); % disp(R2)
                if var(matThetas(1,2,:)) > 1e-10
                    indWsite = nan(J, T1, N); indWtime = nan(J, T1, N); indWj = nan(J, T1, N);
                    for j = 1:J
                        for t = 1:T1
                            for i = 1:N
                                indWsite(j,t,i) = i;  indWtime(j,t,i) = t;  indWj(j,t,i) = j;
                            end
                        end
                    end
                    indWsite = reshape(indWsite, [1,numel(indWsite)]);
                    indWtime = reshape(indWtime, [1,numel(indWtime)]);
                    indWj = reshape(indWj, [1,numel(indWj)]);
                    for j = 1:size(matThetas, 2)
                        subplot(6,6,K0), plot(reshape(matThetas(:,j,:),[nch nsample])');
                        title(strcat('\theta: s =',num2str(indWsite(j)), ',t =', num2str(indWtime(j)), ',j =',num2str(indWj(j)) ))
                        if simu == 1
                            h = refline(0, inittheta(j)); set(h, 'LineWidth',2,'Color','k')
                        end
                        axis tight; %axis off;
                        K0 = K0 +1;
                        if K0 > nbin
                            print('-painters', '-dpsc2', '-r600', '-append', filename)
                            K0 = 1;
                            for ii = 1:nbin
                                subplot(6,6,ii); plot(nan)
                            end
                        end
                    end
                end
            end % monitor theta
            
            % investigating latent variables S
            for monitor_S = 1:monitor_Ses % turn 1:1 to activate
                R3 = psrf(matSes); % disp(R2)
                if var(matSes(1,2,:)) > 1e-10
                    indWsite = nan(T1, N); indWtime = nan(T1, N);
                    for t = 1:T1
                        for i = 1:N
                            indWsite(t,i) = i;  indWtime(t,i) = t;
                        end
                    end
                    indWsite = reshape(indWsite, [1,numel(indWsite)]);
                    indWtime = reshape(indWtime, [1,numel(indWtime)]);
                    for j = 1:size(matSes, 2)
                        subplot(6,6,K0), plot(reshape(matSes(:,j,:), [nch nsample])');
                        title(strcat('S: s =',num2str(indWsite(j)), ',t =', num2str(indWtime(j)) ))
                        if simu == 1
                            h = refline(0, initS(j)); set(h, 'LineWidth',2,'Color','k')
                        end
                        axis tight; %axis off;
                        K0 = K0 + 1;
                        if K0 > nbin
                            print('-painters', '-dpsc2', '-r600', '-append', filename)
                            K0 = 1;
                            for ii = 1:nbin
                                subplot(6,6,ii); plot(nan)
                            end
                        end
                    end
                end
            end
            
            %             % check the posterior predictive distribution
            %             % pay attention to the indices that extracts the components
            %             betahat = reshape(mat(1:Npv,1), [pv,N]);
            %             beta_mean = mat(Npv + (1:p), 1);
            %             if fixdelta0 == 1
            %                 deltahat = 0;
            %                 sigma2hat = mat(npara1, 1);
            %             else
            %                 deltahat = mat(npara1, 1);
            %                 sigma2hat = mat(npara1-1, 1);
            %             end
            %
            %             Yhat = zeros(T,N);
            %             for i = 1:N
            %                 Yhat(:,i) = X(indS==i,varyind)*betahat(:, i) + X(indS==i,comind)*beta_mean(comind) + Smean(2:end,i);
            %             end
            %             Y = reshape(Y,[1,numel(Y)]); Yhat = reshape(Yhat,[1,numel(Yhat)]);
            %             RMSE = sqrt(mean((Y - Yhat).^2));
            %             RSR = sqrt(sum((Y - Yhat).^2))/sqrt(sum((Y - mean(Y)).^2));
            %             NSE = 1 - sum((Y - Yhat).^2)/sum((Y - mean(Y)).^2);
            %             %disp([RMSE, RSR, NSE])
            %
            %             % compute the log-deviance given the "estimated" paramter
            %             tmp = reshape(Y - Yhat, [T, N]);
            %             % for ordinary AR(1) model only
            %             tmp = [tmp(1,:); tmp(2:T,:) - deltahat*tmp(1:(T-1),:)];
            %             tmp = tmp/sqrt(sigma2hat);
            %             L1 = sum(sum( log(normpdf(tmp)) ));
            %             L0 = mean(reshape(L0sAll, [1, numel(L0sAll)]));
            
            %         % compute DIC
            %         D1 = -2*L0; D2 = -2*L1;
            %         DIC = 2*D1 - D2;
            %         pD = D1 - D2;
            %         disp(num2str([D1, pD, DIC]))
            %
            %         % compute DIC for random effects model
            %         D1 = -2*EsAll(1); D2 = -2*EsAll(2); DIC = 2*D1 - D2;
            %         pD = D1 - D2;
            %         disp(num2str([D1, pD, DIC]))
            
            %save(strcat(datapath,'paras',num2str(ncas),'_all.mat'), 'mat','R1','tab1') %,'R2','R3'
            
        end % nvar loop
        
    end
    
end

num2str(Measures,8)

save('paras_pred.mat','Measures','paras','Ys','Y1','namsind')
% disp('done')
end







