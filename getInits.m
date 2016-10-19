function [x] = getInits(x, EV)
% load('paras.mat')
% initvec = paras{1}(:,1)'; %posterior mean of the first model
% save('Inits.mat', 'initvec','namsind')
load('Inits.mat')
N = 67; p = 25; J = 4;
x0 = x;
x0.beta = reshape(initvec(namsind==1), [p,N]);
tmp = initvec(namsind==2)';
x0.beta_mean = tmp(1:p);
x0.beta_sigma2 = tmp(p+(1:p)); 
if EV.common_beta == 1 % all common
    x0.beta_sigma2 = x0.beta_sigma2 * 0;
    x0.beta = repmat(x.beta_mean, [1,N]);
end
x0.sigma2_e = initvec(namsind==3);
x0.delta = initvec(namsind==4);
x0.delta_mean = initvec(namsind==5);
x0.delta_sigma2 = initvec(namsind==6)';
if EV.common_delta == 1 % all common
    x0.delta_sigma2 = x0.delta_sigma2 * 0;
    x0.delta = repmat(x.delta_mean, [1,N]);
end
x0.gamma = reshape(initvec(namsind==7), [J,N]);
tmp = initvec(namsind==8)';
x0.gamma_mean = tmp(1:J);
x0.gamma_sigma2 = tmp(J+(1:J)); 
if EV.common_gamma == 1 % all common
    x0.gamma_sigma2 = x0.gamma_sigma2 * 0;
    x0.gamma = repmat(x.gamma_mean, [1,N]);
end
x0.phi = reshape(initvec(namsind==11), [J,N]);
tmp = initvec(namsind==12)';
x0.phi_mean = tmp(1:J);
x0.phi_sigma2 = tmp(J+(1:J)); 
if EV.common_phi == 1 % all common
    x0.phi_sigma2 = x0.phi_sigma2 * 0;
    x0.phi = repmat(x.phi_mean, [1,N]);
end
tmp = initvec(namsind==9)';
Omega = zeros(J+1);
k = 1; J1 = J+1;
for j = 1:J1
    for i = j:J1
        Omega(i,j) = tmp(k); k = k+1;
        Omega(j,i) = Omega(i,j);
    end
end
x0.Omega = Omega; 
Pre = x0.Omega\eye(J1);
x0.Pre_s = Pre(1,1); % this should = inv(Omega(1,1)-Omega(1,2:J1)*inv(Omega(2:J1,2:J1))*Omega(2:J1,1))
x0.Pre_c = Pre(2:J1, 2:J1); % this should = inv(Omega(2:J1,2:J1)-Omega(2:J1,1)*inv(Omega(1,1))*Omega(1,2:J1))
x0.Ome_s = -Pre(1,2:J1)/Pre(1,1); %should = Omega(1,2:J1)*inv(Omega(2:J1,2:J1))
x0.Ome_c = x0.Omega(2:J1,1)/x0.Omega(1,1);
tmp = initvec(namsind==10)';
Psi = zeros(J);
k = 1;
for j = 1:J
    for i = j:J
        Psi(i,j) = tmp(k); k = k+1;
        Psi(j,i) = Psi(i,j);
    end
end
x0.Psi = Psi; 
x0.invPsi = x0.Psi\eye(J);

x = x0;
end
