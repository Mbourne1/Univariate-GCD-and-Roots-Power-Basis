function [fx_n,gx_n,alpha,theta] = LowRankApprox(fx_n,gx_n,alpha,theta,t)

global SETTINGS

% Perform SNTLN to obtain low rank approximation
switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        fw = GetWithThetas(fx_n,theta);
        gw = GetWithThetas(gx_n,theta);
        
        % Get f(x)+\delta x and g(x) + \delta x for which a low rank
        % approximate Sylvester matrix is obtained.
        
        % Perform STLN - Structured Total Least Norm
        [fw_new,a_gw_new] = STLN(fw,alpha.*gw,t);
        
        fx_n = GetWithoutThetas(fw_new,theta);
               
        % Get g(x) and g(w)
        gw_new = a_gw_new./ alpha;
        gx_n = GetWithoutThetas(gw_new,theta);
        
    case 'STLN Root Specific'
        
        [fw_new,a_gw_new] = STLN_Derivative_Constraint(fx_n,t);

        % Get f(x) = f(x) + \delta f(x)
        fx_n = GetWithoutThetas(fw_new,theta);
        
        % Get g(x) = g(x) + \delta g(x)
        gw_new = a_gw_new./ alpha;
        gx_n = GetWithoutThetas(gw_new,theta);
        
    case 'Standard SNTLN'
        [fx_new,gx_new,alpha_new,theta_new] = SNTLN(fx_n,gx_n,alpha,theta,t);
        
        fx_n = fx_new;
        gx_n = gx_new;
        theta = theta_new;
        alpha = alpha_new;
        
    case 'None'
        
    otherwise
        error('error')
end