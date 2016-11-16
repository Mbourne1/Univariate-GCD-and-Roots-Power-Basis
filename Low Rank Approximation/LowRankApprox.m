function [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = LowRankApprox(fx,gx,alpha,theta,k)
% Get the low rank approximation of the kth Sylvester subresultant matrix
% S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
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
% alpha_lr : 
%
% theta_lr : 
global SETTINGS

% Perform SNTLN to obtain low rank approximation
switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta);
        
        % Get Low Rank Approximation of S_{k}(f,g)
        [fw_lr,a_gw_lr,uw_lr,vw_lr] = STLN(fw,a_gw,k);
        
        
        fx_lr = GetWithoutThetas(fw_lr,theta);
        gx_lr = GetWithoutThetas(a_gw_lr,theta) ./ alpha;
        ux_lr = GetWithoutThetas(uw_lr,theta);
        vx_lr = GetWithoutThetas(vw_lr,theta);
        alpha_lr = alpha;
        theta_lr = theta;
        
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

        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta); 

        [uw,vw] = GetCofactorsCoefficients(fw,a_gw,k);

        ux_lr = GetWithoutThetas(uw,theta);
        vx_lr = GetWithoutThetas(vw,theta);
        
        fx_lr = fx;
        gx_lr = gx;
        
        alpha_lr = alpha;
        theta_lr = theta;
        
    otherwise
        error('error')
end