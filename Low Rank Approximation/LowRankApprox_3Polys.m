function [fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr] ...
    = LowRankApprox_3Polys(fx, gx, hx, alpha, theta, k)
% Get the low rank approximation of the kth Sylvester subresultant matrix
% S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% hx : Coefficients of the polynomial h(x)
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
%
% k : Index of Sylvester subresultant matrix S_{k}(f,g)
%
% % Outputs
%
% fx_lr :
%
% gx_lr :
%
% hx_lr :
%
% ux_lr :
%
% vx_lr :
%
% wx_lr :
%
% alpha_lr :
%
% theta_lr :
global SETTINGS

% Perform SNTLN to obtain low rank approximation
switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        error('Error : Code Branch Not Developed');
        
    case 'STLN Root Specific'
        
        error('Error : Code Branch Not Developed');
        
        
    case 'Standard SNTLN'
        
        error('Error : Code Branch Not Developed');
        
        
    case 'None'
        
        % Get f(\omega) and alpha.*g(\omega) and w(\omega)
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta);
        hw = GetWithThetas(hx,theta);
        
        % Get u(\omega) and v(\omega)
        [uw,vw,ww] = GetCofactorsCoefficients_3Polys(fw,a_gw,hw,k);
        
        % Get u(x) and v(x)
        ux_lr = GetWithoutThetas(uw,theta);
        vx_lr = GetWithoutThetas(vw,theta);
        wx_lr = GetWithoutThetas(ww,theta);
        
        % Get f(x) and g(x) output, which are unchanged from input
        fx_lr = fx;
        gx_lr = gx;
        hx_lr = hx;
        
        % Get \alpha and \theta output, which are unchanged from input
        alpha_lr = alpha;
        theta_lr = theta;
        
    otherwise
        
        error('error')
        
end


end