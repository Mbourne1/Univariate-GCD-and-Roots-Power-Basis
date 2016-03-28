function [fx_n,gx_n,dx, ux, vx, alpha, theta, t , lambda,mu] = o1(fx,gx)
% o1(fx,gx)
%
% Given two polynomials f(x) and g(x) return the GCD d(x)
% 
% Inputs.
%
% fx : Coefficients of polynomial f(x).
%
% gx : Coefficients of polynomial g(x).


% Global variables
global BOOL_PREPROC
global PLOT_GRAPHS
global LOW_RANK_APPROXIMATION_METHOD


% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx) ;


%% Preprocess the polynomials f(x) and g(x)
switch BOOL_PREPROC
    case 'y'
        % Apply Preprocessing
        
        % Get geometric mean
        lambda  = geomean(abs(fx(fx~=0)));
        mu      = geomean(abs(gx(gx~=0)));
        
        % Divide f(x) and g(x) by geometric mean
        fx_n = fx./ lambda;
        gx_n = gx./ mu;
        
        % Get opitmal values of alpha and theta
        [alpha, theta] = GetOptimalAlphaAndTheta(fx_n,gx_n);
        
        % Obtain f(w) and g(w) from f(x) and g(x)]
        fw = GetWithThetas(fx,theta);
        gw = GetWithThetas(gx,theta);
        
        % Plot the unprocessed and preprocessed coefficients of
        % f(x), f(w), g(x) and g(w).
        PlotCoefficients(fx,fw,'f');
        PlotCoefficients(gx,gw,'g');
        
        
        
    case 'n'
        % Dont apply preprocessing
        fx_n = fx;
        gx_n = gx;
        fw = fx;
        gw = gx;
        
        lambda = 1;
        mu = 1;
        
        theta = 1;
        alpha = 1;
        
    otherwise
        error('bool_preproc either y or n');
end

%%
% Get the degree t of the GCD of f(w) and g(w)
t = GetGCDDegree(fw,alpha.*gw);

% Given the degree t, get the optimal column for removal from S_{t}(f,g)
St_preproc  = BuildT(fw,alpha.*gw,t);

% Get the minimum distance between any one of the columns of S_{k}(f,g) and the
% remaining columns of S_{k}(f,g)
[~,opt_col] = GetMinDistance(St_preproc);

%%
% Perform SNTLN to obtain low rank approximation
switch LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        % Get f(x)+\delta x and g(x) + \delta x for which a low rank
        % approximate Sylvester matrix is obtained.
        
        % Perform STLN - Structured Total Least Norm
        [fw_new,a_gw_new] = STLN(fw,alpha.*gw,t);
        
        fx_n = GetWithoutThetas(fw_new,theta);
               
        
        % Get g(x) and g(w)
        gw_new = a_gw_new./ alpha;
        gx_n = GetWithoutThetas(gw_new,theta);
        
    case 'STLN Root Specific'
        
        [fw_new,a_gw_new] = STLN_Derivative_Constraint(fx,t)
        
        % Perform STLN - Structured Total Least Norm
        [fw_new,a_gw_new] = STLN(fw,alpha.*gw,t);
        fx_n = GetWithoutThetas(fw_new,theta);
        
        % Get g(x) and g(w)
        gw_new = a_gw_new./ alpha;
        gx_n = GetWithoutThetas(gw_new,theta);
        
    case 'Standard SNTLN'
        [fx_new,gx_new,alpha_new,theta_new] = SNTLN(fx,gx,t,opt_col,lambda,mu,alpha,theta);
        
        fx_n = fx_new;
        gx_n = gx_new;
        theta = theta_new;
        alpha = alpha_new;
        
        fw_new = GetWithThetas(fx_new,theta);
        gw_new = GetWithThetas(gx_new,theta);
        
        
    case 'None'
        
        fw_new = fw;
        gw_new = gw;
    otherwise
        error('error')
end

St_LowRankApprox = BuildT(fw_new,alpha.*gw_new,t);


% Get the Sylvester matrix of the unprocessed
St_Unproc = BuildT(fx,gx,t);


% Get the quotient polynomials u(x) and v(x)
[ux,vx] = GetQuotients(fx_n,gx_n,t,alpha,theta,opt_col);

% Get the GCD d(x)
dx = GetGCD(ux,vx,fx_n,gx_n,t,alpha,theta);

%%
% Get the singular values for
% 1 : The Sylvester Subresultant for unprocessed f(x) and g(x)
% 2 : The Sylvester Subresultant for preprocessed f(\omega) g(\omega)
% 3 : The Sylvester Subresultant which is a low rank approximation of 2.

vSingularValues_unproc = svd(St_Unproc);
vSingularValues_unproc = Normalise(vSingularValues_unproc);

vSingularValues_preproc = svd(St_preproc);
vSingularValues_preproc = Normalise(vSingularValues_preproc);

vSingularValues_lowRank = svd(St_LowRankApprox);
vSingularValues_lowRank = Normalise(vSingularValues_lowRank);


%%
% Plot Singular values of unproc, preproc, lowrank approx
switch PLOT_GRAPHS
    case 'y'
        figure('name','Singular Values')
        hold on
        plot(log10(vSingularValues_unproc),'-s','DisplayName','Unprocessed')
        plot(log10(vSingularValues_preproc),'-o','DisplayName','Preprocessed')
        plot(log10(vSingularValues_lowRank),'-*','DisplayName','Low Rank Approx')
        legend(gca,'show');
        hold off
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
end

end




