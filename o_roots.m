function [] = o_roots(ex_num,el,bool_preproc,low_rank_approx_method)
% Calculate the roots of a polynomial f(x) by a series of methods.
%
%   Inputs.
%
%
%   ex_num : Example Number
%
%   el : Lower Noise Level

%%
SetGlobals()

%%
% Get polynomial f(x)

EXAMPLE_TYPE = 'FromRoots';

switch EXAMPLE_TYPE 
    case 'FromRoots'
        % Get a set of roots and multiplicities for f(x)
        fx_root_mult_pairs = Examples_Roots(ex_num);
        
        % Get the coefficients of the polynomial f(x)
        fx = GetCoefficients(fx_root_mult_pairs);

        % Print the roots and coefficients of f(x)
        PrintFactorization(fx_root_mult_pairs,'f')
        PrintCoefficientsBivariate(fx,'f')
        
    case 'FromCoefficients'
        switch ex_num
            case '1'
                fx = [1; 1; -2; -2; 1; 1];
            case '2'
                fx = [-1; -1; 1; 1];
            case '3'
                fx = [-5; 1; -5; 1];
        end
    otherwise 
        error('Example polynomial is either from roots or from coefficients')
end


% Add Noise to f(x)
fx = Noise(fx,el);

%%
% Get roots by my method
[computed_root_mult_array] = o_roots_mymethod(fx);

% Analysis.
% Get Distance between polynomial whose roots have been found, and the
% exact input polynomial
fx_computed = GetCoefficients(computed_root_mult_array);
fprintf('Distance between f_exact and f_computed : \n')
display(norm(fx - fx_computed) ./ norm(fx));

% Get roots by matlab method
roots_calc = [roots(flipud(fx)) ones(length(roots(flipud(fx))),1)];
% Print the calculated roots and the corresponding multiplicities.
PrintRoots(roots_calc,'MATLAB METHOD');



% Get roots by zheng method
o_roots_multroot(flipud(fx));

end

