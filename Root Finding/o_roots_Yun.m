function [root_mult_array] = o_roots_Yun(fx)
% Given the polynomial f(x) compute its roots by Square free factorization.
% This algorithm is referred to as Yuns Algorithm in
% https://en.wikipedia.org/wiki/Square-free_polynomial
%
%
% % Inputs
%
% fx : Coefficients of polynomial f(x) as a column vector where first
% coefficient has lowest power [a_{0} ; ... ; a_{m}]^{T}
%
% % Outputs
%
% root_mult_array : A matrix of roots and multiplicities, where the first
% column contains all of the roots of f(x) and the second contains the
% corresponding multiplicities.

global SETTINGS
SETTINGS.PLOT_GRAPHS = 'y';

% Set Iteration number
ite = 1;

fx;
f_dash = Differentiate(fx);

m = GetDegree(fx);
n = GetDegree(f_dash);

deg_limits = [1,min(m,n)];

% Perform GCD computation.

[fx_n,gx_n,u, C{ite}, w{ite}, alpha, theta, t , lambda,mu] ...
    = o_gcd_mymethod(fx,Differentiate(fx),deg_limits);
d{ite} = w{ite}-Differentiate(C{ite});
LineBreakMedium();


while (GetDegree(C{ite}) > 0 )
    
    % Get the degree of polynomial f(x)
    m = GetDegree(C{ite});
    % Get the degree of polynomial g(x)
    n = GetDegree(d{ite});
    
    limits = [1,min(m,n)];
    
    % Get GCD
    [fx,gx,h{ite+1},C{ite+1},w{ite+1},~,~,~,~,~] =...
        o_gcd_mymethod(C{ite},d{ite},limits);
    
   
    
    d{ite+1} = w{ite} - Differentiate(C{ite});
    
    ite = ite+1;
    
    LineBreakMedium();
    
end

SETTINGS.PLOT_GRAPHS = 'n';

root_mult_array = [];

for i = 1:1:length(h)

    try
        
    factor = h{i};
    % Divide by x coefficient
    factor = factor./factor(2);
    % Get root
    root = -1.*factor(1)
    
    catch
        
    end
end

PrintRoots(root_mult_array,'YUNS METHOD');
end