function [ux,vx] = GetQuotients(fx_n,gx_n,t,alpha,theta,opt_col)
% Get Quotient polynomials u(x) and v(x) such that
% f(x)/u(x) = g(x)/v(x) = d(x)
%
%   Inputs.
%
%   fx_n  : Coefficients of polynomial f(x)
%
%   gx_n  : Coefficients of polynomial g(x)
%
%   t     : Degree of GCD d(x)
%
%   alpha : alpha
%
%   theta : theta

% Get degree of polynomial f(x).
m = GetDegree(fx_n);

% Get degree of polynomial g(x).
n = GetDegree(gx_n);

% Get preprocessed form f(\theta,\omega) from f(x) 
fw = GetWithThetas(fx_n,theta);

% Get preprocessed form g(\theta,\omega) from g(x) 
gw = GetWithThetas(gx_n,theta);

% Build the t-th subresultant S_{t}(f,g)
St = BuildT(fw,alpha.*gw,t);

% Get the matrix A_{t}(f,g), S_{t} with the optimal column removed
At = St;
At(:,opt_col) = [];

% Get the column c_{t} removed from S_{t}
ct = St(:,opt_col);

% get the vector x_ls
x_ls = SolveAx_b(At,ct);

% insert 1 into x_ls in the position corresponding to the col
x = [x_ls(1:opt_col-1); -1 ; x_ls(opt_col:end)];

% Split x into v(\omega) and u(\omega)
vw =  x(1:n-t+1);
uw = -x(n-t+2:end);

% Get u(x) from u(w)
ux = GetWithoutThetas(uw,theta);
% Get v(x) from v(w)
vx = GetWithoutThetas(vw,theta);

