function [x] = updateMatrix(x)
global EV N T J J1 C
for myupdateMat = 1:EV.updateMats
    
    A = zeros(J+1); B = zeros(J);
    for i = 1:N
        for t = 2:(T+1)
            w = [x.S(t,i) - x.delta(i)*x.S(t-1,i) - C(:,t-1,i)'*x.gamma(:,i);  C(:,t-1,i)-x.theta(:,t,i)];
            ksi = x.theta(:,t,i) - x.phi(:,i).*x.theta(:,t-1,i);
            if EV.endogeneity == 1
                A = A + w*w';
            else
                A(1,1) = A(1,1) + w(1)^2; A(2:end,2:end) = A(2:end,2:end) + w(2:end)*w(2:end)';
            end
            B = B + ksi*ksi';
        end
    end
    if EV.endogeneity == 1
        x.Omega = iwishrnd(EV.A_Omega+A, EV.nu_Omega);
    else % no endogeneity, Omega(1,2:J1) is fixed to 0
        x.Omega = zeros(J1);
        x.Omega(1,1) = 1./gamrnd( EV.alpha_w,  1/( EV.invbeta_w + 0.5*A(1,1) ) );
        x.Omega(2:end,2:end) = iwishrnd(EV.A_Omega(2:end,2:end)+A(2:end,2:end),  EV.nu_Omega - 1 );
    end
    
    Pre = x.Omega\eye(J1);
    x.Pre_s = Pre(1,1); % this should = inv(x.Omega(1,1)-x.Omega(1,2:J1)*inv(x.Omega(2:J1,2:J1))*x.Omega(2:J1,1))
    x.Pre_c = Pre(2:J1, 2:J1); % this should = inv(x.Omega(2:J1,2:J1)-x.Omega(2:J1,1)*inv(x.Omega(1,1))*x.Omega(1,2:J1))
    x.Ome_s = -Pre(1,2:J1)/Pre(1,1); %should = x.Omega(1,2:J1)*inv(x.Omega(2:J1,2:J1))
    x.Ome_c = x.Omega(2:J1,1)/x.Omega(1,1);
    x.Psi = iwishrnd(EV.A_Psi+B, EV.nu_Psi);
    x.invPsi = x.Psi\eye(J);
    
    % disp(x.Omega)
    
end
end