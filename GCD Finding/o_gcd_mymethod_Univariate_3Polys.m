function [fx_o, gx_o, hx_o, dx_o, ux_o, vx_o, wx_o, alpha_o, theta_o, t , GM_fx, GM_gx, GM_hx] =...
    o_gcd_mymethod_Univariate_3Polys(fx, gx, hx, deg_limits)
% o_gcd_mymethod(fx,gx,deg_limits)
%
% Given two polynomials f(x) and g(x) return the GCD d(x)
%
% % Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% hx : (Vector) Coefficients of polynomial h(x)
%
% deg_limits : (Int Int) Specifiy upper and lower bound of the degree of the GCD
% d(x), typically set when o_gcd_mymethod() is called from a root solving
% problem, where deg_limits have been predetermined.
%
% % Outputs
% 
% fx_o : (Vector) 
%
% gx_o : (Vector)
%
% dx_o : (Vector)
%
% ux_o : (Vector)
% 
% vx_o : (Vector)
%
% alpha_o : (Float)
%
% theta_o : (Float)



% Preprocess the polynomials f(x) and g(x)
[GM_fx, GM_gx, GM_hx, alpha, beta, theta] = Preprocess_3Polys(fx, gx, hx);

% Get f(x), g(x) and h(x) normalised by mean
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;
hx_n = hx./ GM_hx;

% Get f(\omega), g(\omega) and h(\omega)
fw = GetWithThetas(fx_n, theta);
gw = GetWithThetas(gx_n, theta);
hw = GetWithThetas(hx_n, theta);

% Get the degree of the GCD
t = GetGCDDegree_3Polys(fw, alpha.*gw, beta.*hw, deg_limits);
LineBreakLarge();

% Print the degree of the GCD
fprintf([mfilename ' : ' sprintf('Degree of GCD : % i \n',t)]);



if (t == 0)
    dx_o = 1;
    ux_o = fx;
    vx_o = gx;
    return
end


% %
% %
% Get the low rank approximation and refined values for fx, gx, alpha, and
% theta. Note that vectors of polynomial coefficients are suffixed with 
% 'lr' (lr = low rank) and these are the coefficients after low rank
% approximation is computed.

[fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr] = ...
    LowRankApprox_3Polys(fx_n, gx_n, hx_n, alpha, theta, t);


 [ux_lra, vx_lra, wx_lra, fx_lra, gx_lra, hx_lra, dx_lra, alpha_lra, theta_lra] = ...
     APF_3Polys(ux_lr, vx_lr,wx_lr, fx_lr, gx_lr, hx_lr, alpha_lr, theta_lr, t);

% Get f(x), g(x) and h(x)
fx_o = fx_lra;
gx_o = gx_lra;
hx_o = hx_lra;

% Get u(x), v(x) and w(x)
ux_o = ux_lra;
vx_o = vx_lra;
wx_o = wx_lra;

% Get d(x)
dx_o = dx_lra;

% Get \alpha and \theta
alpha_o = alpha_lra;
theta_o = theta_lra;
end




