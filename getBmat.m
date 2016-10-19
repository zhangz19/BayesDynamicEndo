function [out] = getBmat(phi, T, varargin)
% get information of AR(1)-type correlation matrix
mat_ = any(strcmp(varargin, 'mat'));
inv_ = any(strcmp(varargin, 'inv'));
chol_ = any(strcmp(varargin, 'chol'));
det_ = any(strcmp(varargin, 'det'));
time_ = any(strcmp(varargin, 'time'));

t0 = cputime;

if mat_
    out.mat = zeros(T);
    for t1 = 1:T
        for t2 = 1:T
            out.mat(t1,t2) = phi^(abs(t1-t2));
        end
    end
end

if inv_
    out.inv = zeros(T); h = 1/(1-phi^2); hphi = h*phi; hp2 = h*(1+phi^2);
    out.inv(1,1) = h; out.inv(1,2) = -hphi;
    out.inv(T,T) = h; out.inv(T,T-1) = -hphi;
    for t1 = 2:(T-1)
        out.inv(t1,t1) = hp2; 
        out.inv(t1,t1-1) = -hphi; out.inv(t1,t1+1) = -hphi;
    end
end

if chol_
    out.chol = zeros(T); out.chol(:,1) = phi.^(0:(T-1)); h = sqrt(1-phi^2);
    for t1 = 1:(T-1)
        tmp = phi^(t1-1)*h;
        for t2 = 1:(T-t1)
            out.chol(t2+t1,t2+1) = tmp;
        end
    end
end

if det_
    out.det = (1-phi^2)^(T-1);
end

if time_
    disp('The cputime is:'); disp(cputime - t0)
end