function [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = LowRankApprox_2Polys(fx, gx, alpha, theta, k)
% Get the low rank approximation of the kth Sylvester subresultant matrix
% S_{k}(f,g)
%
% % Inputs
%
% [fx, gx] : Coefficients of the polynomial f(x) and g(x)
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
%
% k : Index of Sylvester subresultant matrix S_{k}(f,g)
%
% % Outputs
%
% [fx_lr, gx_lr] : Coefficients of f(x) and g(x)
%
% [ux_lr, vx_lr] : Coefficients of u(x) and v(x)
%
% alpha_lr :
%
% theta_lr :
global SETTINGS

% Perform SNTLN to obtain low rank approximation
switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        % Get f(\omega) and \alpha.*g(\omega)
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta);
        
        % Get Low Rank Approximation of S_{k}(f,g)
        [fw_lr, a_gw_lr, uw_lr, vw_lr] = STLN(fw, a_gw, k);
        
        % Get f(x) and g(x) output
        fx_lr = GetWithoutThetas(fw_lr,theta);
        gx_lr = GetWithoutThetas(a_gw_lr,theta) ./ alpha;
        
        % Get u(x) and v(x) output
        ux_lr = GetWithoutThetas(uw_lr,theta);
        vx_lr = GetWithoutThetas(vw_lr,theta);
        
        % Get \alpha and \theta output
        alpha_lr = alpha;
        theta_lr = theta;
        
        
        if(SETTINGS.PLOT_GRAPHS)
            
                S1 = BuildT(fx,gx,k);
                S2 = BuildT(fw,a_gw,k);
                S3 = BuildT(fx_lr,gx_lr,k);
                S4 = BuildT(fw_lr,a_gw_lr,k);
                
                vSingularValues1 = svd(S1);
                vSingularValues2 = svd(S2);
                vSingularValues3 = svd(S3);
                vSingularValues4 = svd(S4);
                
                figure()
                plot(log10(vSingularValues1),'-s','DisplayName','f(x) g(x)');
                hold on
                plot(log10(vSingularValues2),'-s','DisplayName','f(\omega) g(\omega)');
                plot(log10(vSingularValues3),'-s','DisplayName','f(x)_lr g(x)_lr');
                plot(log10(vSingularValues4),'-s','DisplayName','f(\omega)_lr g(\omega)_lr');
                hold off
            
        end
        
    case 'STLN Root Specific'
        
        error('err - Not tested this branch');
        % Get f(\omega)
        fw = GetWithThetas(fx,theta);
        
        [fw_lr,a_gw_lr,uw_lr,vw_lr] = STLN_Derivative_Constraint(fw,k);
        
        % Get f(x) = f(x) + \delta f(x)
        fx_lr = GetWithoutThetas(fw_lr,theta);
        gx_lr = GetWithoutThetas(a_gw_lr,theta)./alpha;
        
        ux_lr = GetWithoutThetas(uw_lr,theta);
        vx_lr = GetWithoutThetas(vw_lr,theta);
        
        theta_lr = theta;
        alpha_lr = alpha;
        
    case 'Standard SNTLN'
        
        [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = SNTLN(fx,gx,alpha,theta,k);
        
        
        
        
    case 'None'
        
        % Get f(\omega) and alpha.*g(\omega)
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta);
        
        % Get u(\omega) and v(\omega)
        [uw,vw] = GetCofactorsCoefficients_2Polys(fw,a_gw,k);
        
        % Get u(x) and v(x)
        ux_lr = GetWithoutThetas(uw,theta);
        vx_lr = GetWithoutThetas(vw,theta);
        
        % Get f(x) and g(x) output, unchanged from input
        fx_lr = fx;
        gx_lr = gx;
        
        % Get \alpha and \theta output, unchanged from input
        alpha_lr = alpha;
        theta_lr = theta;
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
    otherwise
        
        error('error')
        
end


end