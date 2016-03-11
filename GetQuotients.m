function [ux,vx] = GetQuotients(fx_n,gx_n,t,alpha,theta,opt_col)
% Get Quotient polynomials u(x) and v(x) such that
% f(x)/u(x) = g(x)/v(x) = d(x)
%
%   Inputs.
%
%   fx_n  :
%
%   gx_n  :
%
%   t     :
%
%   alpha :
%
%   theta :

% Get degree of polynomial f(x).
[rows_f,~] = size(fx_n);
m = rows_f - 1;

% Get degree of polynomial g(x).
[rows_g,~] = size(gx_n);
n = rows_g - 1;

% Get preprocessed form f(\theta,\omega) g(\theta,\omega)
fw = fx_n.*(theta.^(0:1:m)');
gw = gx_n.*(theta.^(0:1:n)');

% Build the t-th subresultant S_{t}(f,g)
St = BuildSubresultant(fw,gw,t,alpha);


% Get the matrix A_{t}(f,g), S_{t} with the optimal column removed
At = St;
At(:,opt_col) = [];

% Get the column c_{t} removed from S_{t}
ct = St(:,opt_col);

% get the vector x_ls
x_ls = SolveAx_b(At,ct);

% Get associated residual
res = norm(ct - (At*x_ls));

% insert 1 into x_ls in the position corresponding to the col
x = [x_ls(1:opt_col-1); -1 ; x_ls(opt_col:end)];

% Split x into v(\omega) and u(\omega)
vw =  x(1:n-t+1);
uw = -x(n-t+2:end);


ux = uw./(theta.^((0:1:m-t)'));
vx = vw./(theta.^((0:1:n-t)'));
