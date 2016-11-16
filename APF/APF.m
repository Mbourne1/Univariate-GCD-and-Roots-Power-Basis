function [ux_lra, vx_lra, fx_lra, gx_lra, dx_lra, alpha_lra, theta_lra] ...
    = APF(ux, vx, fx, gx, alpha, theta, k)
% Get the coefficients of the GCD d(x), of f(x) and g(x), given u(x) and v(x).
% This function has two branches. d(x) can be computed either by utilising
% [u(x), v(x),f(x) and g(x)], or just [u(x) and f(x)]
%
% % Inputs.
%
% ux : Coefficients of input polynomial u(x)
%
% vx : Coefficients of input polynomial v(x)
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
%
% k : Degree of the GCD d(x)

global SETTINGS
switch SETTINGS.APF_METHOD
    case 'None'
        
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta);
        uw = GetWithThetas(ux,theta);
        vw = GetWithThetas(vx,theta);
        
        [dw] = GetGCDCoefficients(uw, vw, fw, a_gw,k);
        
        dx_lra = GetWithoutThetas(dw,theta);
        ux_lra = ux;
        vx_lra = vx;
        fx_lra = fx;
        gx_lra = gx;
        alpha_lra = alpha;
        theta_lra = theta;
        
    case 'Standard APF Nonlinear'
        
        error([mfilename ' : ' 'Code not completed']);
        
    case 'Standard APF Linear'
        
        error([mfilename ' : ' 'Code not completed']);
        
    otherwise
        error('err')
end



end
