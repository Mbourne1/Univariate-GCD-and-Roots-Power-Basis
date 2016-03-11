function [dx] = GetGCD(ux,vx,fx,gx,t,alpha,theta)
% Get the GCD of f(x) and g(x), given u(x) and v(x).
%
%   Inputs.
%
%   ux      :
%
%   vx      : 
%
%   fx      :
%
%   gx      :
%
%   t       : 
%
%   alpha   :
%
%   theta   :

% Get degree of polynomial f(x)
[rows_f,~] = size(fx);
m = rows_f -1;

% Get degree of polynomial g(x)
[rows_g,~] = size(gx);
n = rows_g -1;

% Get u(w) and v(w)
uw = ux .* (theta.^(0:1:m-t)');
vw = vx .* (theta.^(0:1:n-t)');

% Get f(w) and g(w)
fw = fx .* theta.^(0:1:m)';
gw = gx .* theta.^(0:1:n)';

% Get the GCD d(x) by the APF
C_u = BuildC1(uw,t);
C_v = BuildC1(vw,t);

C1 = ...
    [
    C_u;
    C_v;
    ];

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = [fw;alpha.*gw];

% Get the vector d.
d_ls = SolveAx_b(C1,rhs_vec)

% Calculate d(\omega)
dw = d_ls;

% Get d(x)
dx = dw./(theta.^(0:1:t))';
