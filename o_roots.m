function [] = o_roots(ex_num,el,bool_preproc,low_rank_approx_method)
% Calculate the roots of a polynomial f(x) by a series of methods.
%
% Inputs.
%
% ex_num : Example Number
%
% el : Lower Noise Level
% 
% bool_preproc :
% 
% low_rank_approx_method : 

% %
% Set global variables.
SetGlobalVariables(bool_preproc,low_rank_approx_method)

%%
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

