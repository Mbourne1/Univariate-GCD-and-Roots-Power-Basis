function [ux,vx] = GetCofactorsCoefficients(fx,gx,t,opt_col)
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

% Get degree of polynomial g(x).
n = GetDegree(gx);

% Build the t-th subresultant S_{t}(f,g)
St = BuildT(fx,gx,t);

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
vx =  x(1:n-t+1);
ux = -x(n-t+2:end);


end

