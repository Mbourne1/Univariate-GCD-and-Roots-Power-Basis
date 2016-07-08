function [root_mult_array] = o_roots_Musser(fx)
% Given the polynomial f(x) compute its roots by Square free factorization.
% This algorithm is referred to as Musser's Algorithm in 
% http://www.cs.berkeley.edu/~fateman/282/F%20Wright%20notes/week6.pdf 
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


% Set Iteration number
ite = 1;

% 
f{1} = fx;

% Get the derivative of f(x)
gx{1} = Differentiate(f{1});

% Get the degree of f(x)
m = GetDegree(f{1});

n = GetDegree(gx{1});

% Set upper and lower limits of GCD(f,g) # Since number of distinct roots
% is unknown, upper and lower limits are unknown.
lower_lim = 1;
upper_lim = min(m,n);
deg_limits = [lower_lim, upper_lim];

% Perform GCD computation.
[fx_n,gx_n,dx, ux_o, vx_o, alpha, theta, t , lambda,mu] ...
    = o_gcd_mymethod(f{1},gx{ite},deg_limits);


LineBreakMedium();
g{ite} = dx;

h{ite} = Deconvolve(f{ite},g{ite});

while (GetDegree(h{ite}) > 0 )
    
    % Get the degree of polynomial f(x)
    m = GetDegree(g{ite});
    % Get the degree of polynomial g(x)
    n = GetDegree(h{ite});
    
    % Set Limits
    lower_lim = 1;
    upper_lim = min(m,n);
    deg_limits = [lower_lim, upper_lim];
    
    if (GetDegree(h{ite}) ==0 || GetDegree(g{ite}) == 0)
        h{ite+1} = 1;
        vx{ite} = h{ite};
        ux{ite} = g{ite};
    else
        
          [fx_n,gx_n,h{ite+1}, ux{ite}, vx{ite}, ~, ~, ~ , ~,~] ...
              = o_gcd_mymethod(g{ite},h{ite},deg_limits);
        
    end
    
    % The polynomial g can be obtained in two ways, as u(x) from the GCD
    % triple (d(x),u(x),v(x)) or by deconvolution.
   
   
    switch SETTINGS.ROOTS_UX
        case 'From ux'
            g{ite+1} = ux{ite};
        case 'From Deconvolution'
            g{ite+1} = Deconvolve(g{ite},h{ite+1});
        otherwise
            error('err')
    end
    
    ite = ite+1;
    
    LineBreakMedium();
    
end

root_mult_array = [];

for i = 1:1:length(vx)
    
    try
        %fprintf('Roots of multiplicity %i \n',i)
        
        factor = vx{i};
        m = GetDegree(factor);
        if m > 1
            % Get roots by matlab method
            root = roots(flipud(factor));
            
            % Get number of roots
            nRoots = length(root);
            
            % Get the new roots and mutliplicities
            new_roots_mults = [root i.*ones(nRoots,1)];
            
            % Add the new roots to the array of roots
            root_mult_array = [root_mult_array ; new_roots_mults];
            
        else
            % Divide by x coefficient
            factor = factor./factor(2);
            % Get root
            root = -1.*factor(1);
            
            new_root_mult = [root i];
            
            % Add the new root to the array of roots
            root_mult_array = [root_mult_array ; new_root_mult];
            
            
        end

        
    catch
    end
end

PrintRoots(root_mult_array,'MUSSER METHOD');
end