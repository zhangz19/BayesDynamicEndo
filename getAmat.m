function [out] = getAmat(tau, rho, DistMat_r, varargin)
% get information of exponential correlation matrix
% turn tau = 1 to get information for D part only
mat = tau*exp(-rho*DistMat_r);

mat_ = any(strcmp(varargin, 'mat'));
inv_ = any(strcmp(varargin, 'inv'));
chol_ = any(strcmp(varargin, 'chol'));
cholinv_ = any(strcmp(varargin, 'cholinv'));
det_ = any(strcmp(varargin, 'det'));
time_ = any(strcmp(varargin, 'time'));

t0 = cputime;
if inv_ && ~chol_ && ~det_ % inverse only
    out.inv = mat\eye(size(mat, 1)); %best way in MATLAB to inverse matrix
elseif chol_ || det_ || cholinv_
    cholmat = chol(mat,'lower');
    if det_
        out.det = (prod(diag(cholmat)))^2;
    end
    if cholinv_ || inv_
        out.cholinv = cholmat\eye(size(cholmat,1));
        if inv_
            out.inv = cholmat'\out.cholinv;
        end
    end
    if chol_
       out.chol = cholmat;
    end
    if ~cholinv_ && inv_
       out = rmfield(out, 'cholinv');
    end
end
if mat_
    out.mat = mat;
end
if time_
    disp('The cputime is:'); disp(cputime - t0)
end