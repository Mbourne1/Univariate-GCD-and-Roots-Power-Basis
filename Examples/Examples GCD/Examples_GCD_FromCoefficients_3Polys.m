function [fx, gx, hx, dx, ux, vx, wx] = Examples_GCD_FromCoefficients_3Polys(ex_num)
%
% % Inputs
%
% ex_num : Example Number
%
% % Outputs
%
% [fx, gx, hx] : Coefficients of the polynomials f(x), g(x) and h(x)
%
% dx : Coefficients of the common divisor d(x)
%
% [ux, vx, wx] : Coefficients of the polynomials u(x), v(x) and w(x)

addpath(genpath('../Examples'))

% Get the rm (Root - Multiplicity) array for each of the polynomials
[f_rm_array, g_rm_array, h_rm_array, d_rm_array, u_rm_array, v_rm_array, w_rm_array] = ...
    GCD_Examples_Univariate_3Polys(ex_num);


% Get coefficients of the polynomials f(x), g(x) and h(x)
fx = GetCoefficientsFromSymbolicRoots(f_rm_array);
gx = GetCoefficientsFromSymbolicRoots(g_rm_array);
hx = GetCoefficientsFromSymbolicRoots(h_rm_array);

% Get coefficients of the GCD d(x)
dx = GetCoefficientsFromSymbolicRoots(d_rm_array);

% Get coefficients of u(x), v(x) and w(x)
ux = GetCoefficientsFromSymbolicRoots(u_rm_array);
vx = GetCoefficientsFromSymbolicRoots(v_rm_array);
wx = GetCoefficientsFromSymbolicRoots(w_rm_array);

fx_sym = GetSymbolicPolyFromSymbolicRoots(f_rm_array);
gx_sym = GetSymbolicPolyFromSymbolicRoots(g_rm_array);
hx_sym = GetSymbolicPolyFromSymbolicRoots(h_rm_array);

dx_sym = GetSymbolicPolyFromSymbolicRoots(d_rm_array);

ux_sym = GetSymbolicPolyFromSymbolicRoots(u_rm_array);
vx_sym = GetSymbolicPolyFromSymbolicRoots(v_rm_array);
wx_sym = GetSymbolicPolyFromSymbolicRoots(w_rm_array);

display(fx_sym)
display(gx_sym)
display(hx_sym)

display(dx_sym)

display(ux_sym)
display(vx_sym)
display(wx_sym)

end

