function [ux_lra, vx_lra, wx_lra, fx_lra, gx_lra, hx_lra, dx_lra, alpha_lra, theta_lra] ...
    = APF_3Polys(ux, vx, wx, fx, gx, hx, alpha, theta, k)
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
% wx :
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% hx : 
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
%
% k : Degree of the GCD d(x)
%
% % Outputs
%
% [ux_lra vx_lra wx_lra]
%
% [fx_lra gx_lra hx_lra]
%
% dx_lra
%
% alpha_lra
%
% theta_lra

global SETTINGS
switch SETTINGS.APF_METHOD
    case 'None'
        
        % Get f(\omega) and g(\omega) and h(\omega)
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta);
        hw = GetWithThetas(hx,theta);
        
        % Get u(\omega) and v(\omega) and w(\omega)
        uw = GetWithThetas(ux,theta);
        vw = GetWithThetas(vx,theta);
        ww = GetWithThetas(wx,theta);
        
        % Get d(\omega)
        [dw] = GetGCDCoefficients_3Polys(uw, vw, ww, fw, a_gw, hw, k);
        
        % Get d(x) output
        dx_lra = GetWithoutThetas(dw,theta);
        
        % Get u(x) and v(x) output unchanged from input
        ux_lra = ux;
        vx_lra = vx;
        wx_lra = wx;
        
        % Get f(x) and g(x) output, unchanged from input
        fx_lra = fx;
        gx_lra = gx;
        hx_lra = hx;
        
        % Get \alpha and \theta output, unchanged from input
        alpha_lra = alpha;
        theta_lra = theta;
        
        SETTINGS.APF_REQ_ITE = 0;
        
    case 'Standard APF Nonlinear'
        
        error([mfilename ' : ' 'Code not completed']);
        SETTINGS.APF_REQ_ITE = 0;
    case 'Standard APF Linear'
        
        error([mfilename ' : ' 'Code not completed']);
        SETTINGS.APF_REQ_ITE = 0;
    otherwise
        error('err')
end



end
