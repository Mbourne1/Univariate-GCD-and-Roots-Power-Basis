function [fx_n,gx_n,dx, ux, vx, alpha, theta, t , lambda,mu] =...
    o1(fx,gx,nDistinctRoots)
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
global PLOT_GRAPHS


% Preprocess the polynomials f(x) and g(x)
[lambda,mu,alpha,theta] = Preproccess(fx,gx);

% Get f(x) normalised by mean
fx_n = fx./lambda;

% Get g(x) normalised by mean
gx_n = gx./mu;

% Get f(w)
fw = GetWithThetas(fx_n,theta);

% Get g(w)
gw = GetWithThetas(gx_n,theta);

% %
% %
% Get the degree t of the GCD of f(w) and g(w)
if (nDistinctRoots == 1)
    fprintf('Only one distinct root \n')
    fprintf('Degree of GCD is equal to g(x)')
end

% %
% %
% Get the degree of the GCD
t = GetGCDDegree(fw,alpha.*gw,nDistinctRoots);

% Given the degree t, get the optimal column for removal from S_{t}(f,g)
St_preproc  = BuildT(fw,alpha.*gw,t);

% Build the Sylvester subresultant matrix for the unprocessed polynomials
% f(x) and g(x)
S1_preproc = BuildT(fw,alpha.*gw,1);

% Get the minimum distance between any one of the columns of S_{k}(f,g) and the
% remaining columns of S_{k}(f,g)
[~,opt_col] = GetMinDistance(St_preproc);



% %
% %
% Get the low rank approximation and refined values for fx,gx,alpha, and
% theta.

[fx_n,gx_n,alpha,theta] = LowRankApprox(fx_n,gx_n,alpha,theta,t);
fw_new = GetWithThetas(fx_n,theta);
gw_new = GetWithThetas(gx_n,theta);

% Build the Sylvester matrix which is the low rank approximation
S1_LowRankApprox = BuildT(fw_new,alpha.*gw_new,1);

% Build the Sylvester matrix of the unprocessed input polynomials.
S1_Unproc = BuildT(fx,gx,1);

% %
% Get the quotient polynomials u(x) and v(x)
[ux,vx] = GetQuotients(fx_n,gx_n,t,alpha,theta,opt_col);

% %
% Get the GCD d(x)
dx = GetGCD(ux,vx,fx_n,gx_n,t,alpha,theta);

% %
% %
% %
% Get the singular values for
% 1 : The Sylvester Subresultant for unprocessed f(x) and g(x)
% 2 : The Sylvester Subresultant for preprocessed f(\omega) g(\omega)
% 3 : The Sylvester Subresultant which is a low rank approximation of 2.

vSingularValues_unproc = svd(S1_Unproc);
vSingularValues_unproc = Normalise(vSingularValues_unproc);

vSingularValues_preproc = svd(S1_preproc);
vSingularValues_preproc = Normalise(vSingularValues_preproc);

vSingularValues_lowRank = svd(S1_LowRankApprox);
vSingularValues_lowRank = Normalise(vSingularValues_lowRank);


% %
% %
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




