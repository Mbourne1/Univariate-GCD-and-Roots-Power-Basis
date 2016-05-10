function [fx] = Examples_Roots(ex_num)


EXAMPLE_TYPE = 'From Roots';

switch EXAMPLE_TYPE
    case 'From Roots'
        
        % Get a set of roots and multiplicities for f(x)
        fx_root_mult_array = Examples_Roots_FromRoots(ex_num);
        
        % Get the coefficients of the polynomial f(x)
        fx = GetCoefficients(fx_root_mult_array);
        
        % Print the roots and coefficients of f(x)
        PrintFactorization(fx_root_mult_array,'f')
        PrintCoefficientsBivariate(fx,'f')
        
    case 'From Coefficients'
        
        % Get the coefficients of polynomial f(x)
        fx = Examples_Roots_FromCoefficients(ex_num);
        PrintCoefficientsBivariate(fx,'f')
        
    otherwise
        error('Example polynomial is either from roots or from coefficients')
end

end


function root_mult_arr = Examples_Roots_FromRoots(ex_num)
% Given an example number, return a polynomial in terms of its roots and
% their corresponding multiplicities.
%
% Input
%
% Example Number
%
% Outputs.
%
% root_mult_arr : Array of roots and corresponding multiplicities in f(x)

switch ex_num
    case '1'
        root_mult_arr = ...
            [
            0.1 7;
            0.9 12;
            ];
    case '2'
        root_mult_arr = ...
            [
            0.1 1;
            0.2 2;
            0.9 3;
            ];
        
    case '3'
        root_mult_arr = ...
            [
            6   10;
            -3  3;
            -7  1;
            ];
    case '4'
        root_mult_arr = ...
            [
            1   1
            2   2
            3   3
            4.5 1
            ];
    case '5'
        root_mult_arr = ...
            [
            1.05467 1
            2.24587 2
            5.54743 3
            1.75647 2
            ];
    case '6'
        root_mult_arr = ...
            [
            1.05467 1
            2.24587 2
            5.54743 3
            1.75647 2
            2.56478 5
            1.15445 1
            ];
    case 'Custom'
        prompt = 'Enter the degree of Polynomial f(x) :';
        m = input(prompt);
        
        intvl_low = -1;
        intvl_high = +1;
        
        root_mult_arr = BuildRandomPolynomial(m,intvl_low,intvl_high);
        
    case 'Wilkinsons'
        root_mult_arr = ...
            [
            1 1
            2 1
            3 1
            4 1
            5 1
            6 1
            ];
    case 'Zeng'
        root_mult_arr = ...
            [
            1   20
            2   15
            3   10
            4   5
            ];
        
    case 'Zeng2'
        % A Method Computing Multiple Roots of Inexact Polynomials - Zeng
        % Dayton 
        root_mult_arr = ...
            [
            1   4
            2   3
            3   2
            4   1
            ];
        
    otherwise
        error('err')
end




end


function fx_exact = Examples_Roots_FromCoefficients(ex_num)

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

end