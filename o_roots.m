function [] = o_roots(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_roots(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
%
% Calculate the roots of a polynomial f(x) by a series of methods.
%
% Inputs.
%
% ex_num : Example Number
%
% el : Lower Noise Level
% 
% mean_method : Method of mean used to divide coefficients of Sylvester 
%               Matrix as part of preprocessing.
% 'None' : No mean method is used
% 'Geometric Mean Matlab Method' : Divide by Geometric Mean
%
% bool_alpha_theta : 
%   'y' : Include Preprocessing
%   'n' : Exclude Preprocessing
% 
% low_rank_approx_method : 
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
%   'Root Specific SNTLN'
%
%
% % Example
% o_roots('1',1e-10,'Geometric Mean Matlab Method', 'y','Standard STLN')
% 


% %
% Set global variables.
SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)

global GCD_OR_ROOTS
GCD_OR_ROOTS = 'Roots';

% %
% Get polynomial f(x)

EXAMPLE_TYPE = 'FromRoots';
%EXAMPLE_TYPE = 'FromCoefficients'

switch EXAMPLE_TYPE 
    case 'FromRoots'
        % Get a set of roots and multiplicities for f(x)
        fx_root_mult_pairs = Examples_Roots(ex_num);
        
        % Get the coefficients of the polynomial f(x)
        fx_exact = GetCoefficients(fx_root_mult_pairs);

        % Print the roots and coefficients of f(x)
        PrintFactorization(fx_root_mult_pairs,'f')
        PrintCoefficientsBivariate(fx_exact,'f')
        
    case 'FromCoefficients'
        switch ex_num
            case '1'
                fx_exact = [1; 1; -2; -2; 1; 1];
            case '2'
                fx_exact = [-1; -1; 1; 1];
            case '3'
                fx_exact = [-5; 1; -5; 1];
            case '4'
                fx_exact = [48; -32; 0; 0; 1];
        end
    otherwise 
        error('Example polynomial is either from roots or from coefficients')
end


% Add Noise to f(x)
fx = Noise(fx_exact,el);

%% MY METHOD
% Get roots by my method and compare the computed f(x) with the exact f(x)
[computed_root_mult_array_mymethod] = o_roots_mymethod(fx);


fx_computed_mymethod = GetCoefficients(computed_root_mult_array_mymethod);
fprintf('Distance between f_exact and f_computed mymethod : \n')
display(norm(fx_exact - fx_computed_mymethod) ./ norm(fx_exact));

%% MATLAB
% Get roots by matlab method

roots_calc_matlab = o_roots_matlab(fx);
fx_computed_matlab = GetCoefficients(roots_calc_matlab);
fprintf('Distance between f_exact and f_computed matlab : \n')
display(norm(fx_exact - fx_computed_matlab) ./ norm(fx_exact));


%% MULTROOT
% Get roots by zheng method
computed_root_mult_array_multroot = o_roots_multroot(flipud(fx));
fx_computed_multroot = GetCoefficients(computed_root_mult_array_multroot);
fprintf('Distance between f_exact and f_computed multroot method: \n')
display(norm(fx_exact - fx_computed_multroot) ./ norm(fx_exact));


end

