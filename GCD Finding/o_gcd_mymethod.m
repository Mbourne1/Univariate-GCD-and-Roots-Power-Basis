function [fx_o, gx_o, dx_o, ux_o, vx_o, alpha_o, theta_o, t , lambda,mu] =...
    o_gcd_mymethod(fx, gx, deg_limits)
% o_gcd_mymethod(fx,gx,deg_limits)
%
% Given two polynomials f(x) and g(x) return the GCD d(x)
%
% % Inputs.
%
% fx : Coefficients of polynomial f(x).
%
% gx : Coefficients of polynomial g(x).
%
% deg_limits : Specifiy upper and lower bound of the degree of the GCD
% d(x), typically set when o_gcd_mymethod() is called from a root solving
% problem, where deg_limits have been predetermined.
%
% % Outputs
% 
% fx_o : Coefficients of polynomial f(x)
%
% gx_o : Coefficients of polynomial g(x)
%
% dx_o : Coefficients of polynomial d(x)
%
% ux_o : Coefficients of polynomial u(x)
% 
% vx_o : Coefficients of polynomial v(x)
%
% alpha_o : Optimal value of \alpha
%
% theta_o : Optimal value of \theta



% Preprocess the polynomials f(x) and g(x)
[lambda,mu,alpha,theta] = Preprocess(fx,gx);

% Get f(x) normalised by mean
fx_n = fx./lambda;

% Get g(x) normalised by mean
gx_n = gx./mu;

% Get f(\omega)
fw = GetWithThetas(fx_n,theta);

% Get g(\omega)
gw = GetWithThetas(gx_n,theta);


% Get the degree of the GCD with limits defined
t = GetGCDDegree(fw,alpha.*gw,deg_limits);
LineBreakLarge();

% Print the degree of the GCD
fprintf([mfilename ' : ' sprintf('Degree of GCD : %s \n',int2str(t))]);



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

[fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = LowRankApprox(fx_n,gx_n,alpha,theta,t);


[ux_lra, vx_lra, fx_lra, gx_lra, dx_lra, alpha_lra, theta_lra] = ...
    APF(ux_lr, vx_lr, fx_lr, gx_lr, alpha_lr, theta_lr, t);

% Get f(x) and g(x)
fx_o = fx_lra;
gx_o = gx_lra;

% Get u(x) and v(x)
ux_o = ux_lra;
vx_o = vx_lra;

% Get d(x)
dx_o = dx_lra;

% Get \alpha and \theta
alpha_o = alpha_lra;
theta_o = theta_lra;
end




