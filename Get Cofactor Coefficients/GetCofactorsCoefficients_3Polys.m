function [ux, vx, wx] = GetCofactorsCoefficients_3Polys(fx, gx, hx, k)
% Get Quotient polynomials u(x) and v(x) such that
% f(x)/u(x) = g(x)/v(x) = d(x)
%
% % Inputs.
%
% fx  : Coefficients of polynomial f(x)
%
% gx  : Coefficients of polynomial g(x)
%
% hx  : Coefficients of polynomial h(x)
%
% k     : Degree of common divisor d(x) of degree k.
%
% % Outputs
%
% ux : Coefficients of the polynomial u(x)
%
% vx : Coefficients of the polynomial v(x)
%
% wx : Coefficients of the polynomial w(x)


% Get degree of polynomial g(x).
n = GetDegree(gx);

% Get degree of polynomail h(x)
o = GetDegree(hx);

% Build the t-th subresultant S_{t}(f,g)
Sk = BuildT_3Polys(fx, gx, hx, k);

% Get index of optimal column for removal from S_{k}(f,g)
[~, idx_col] = GetMinDistance(Sk);

% Get the matrix A_{t}(f,g), S_{t} with the optimal column removed
Ak = Sk;
Ak(:,idx_col) = [];

% Get the column c_{t} removed from S_{t}
ck = Sk(:,idx_col);

% Get the vector x_ls
x_ls = SolveAx_b(Ak,ck);

% Insert 1 into x_ls in the position corresponding to the col
x = [x_ls(1:idx_col-1); -1 ; x_ls(idx_col:end)];


% Get the number of coefficients in v(x)
nCoeffs_vx = n-k+1;

% Get the number of coefficients in w(x)
nCoeffs_wx = o-k+1;

% Get coefficients of v(x), w(x) and u(x)
vx =  x(1:nCoeffs_vx);
wx =  x(nCoeffs_vx+1 : nCoeffs_vx+nCoeffs_wx)
ux = -x(nCoeffs_vx+nCoeffs_wx + 1:end);


end

