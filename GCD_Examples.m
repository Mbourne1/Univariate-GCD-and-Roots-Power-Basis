
function [root_mult_arr_fx,root_mult_arr_gx,roots_dx,roots_ux,roots_vx,rankloss] = GCD_Examples(n)
% Given an example number. 
%
% Inputs.
%
% n : Index of example to be used
%
% Outputs.
%
% a : Roots and multiplicities of polynomial f.
%
% b : Roots and multiplicities of polynomial g.
%
% c : Roots and multiplicities of polynomial d, the GCD of f and g.
%
% u : Roots and multiplicities of quotient polynomial f/u = d.
%
% v : Roots and multiplicities of quotient polynomial g/v = d.
%

switch n
    
    case '1' % From Bini
        root_mult_arr_fx = ...
            [
                -1      1;
                -2     1;
                2      1;
                1      1;
            ];
        
        root_mult_arr_gx = ...
            [
                2      1;
                -0.5   1;
                0.5    1;
            ];
        
    case '2'
        root_mult_arr_fx = [
            -1      1;
            -2     1;
            2      1;
            1      1;
            3      1;
            -3     1;
            ];
        
        root_mult_arr_gx = [
            2      1;
            -0.5   1;
            0.5    1;
            3      1;
            -3     1;
            ];
        
    case '3'
        root_mult_arr_fx = [
            0.10    10
            0.56    4
            0.40    4
            0.79    3
            0.69    2
            ];
        
        root_mult_arr_gx = [
            0.10    10
            0.56    4
            0.69    2
            ];
        
        
        
    case '4'
        root_mult_arr_fx = [
            0.2 2
            0.5 1
            0.7 1
            0.9 1
            1.1 1
            ];
        root_mult_arr_gx = ...
            [
            0.5 1
            0.1 1
            0.9 1
            0.3 1
            ];
        
        
    case '5'
        root_mult_arr_fx = ...
            [
            0.2 1
            0.4 1
            0.6 1
            0.8 1
            ];
        
        root_mult_arr_gx = ...
            [
            0.9 1
            0.2 1
            0.3 1
            ];
        
    case '6'
        root_mult_arr_fx = [
            0.56 20
            0.75 3
            0.82 3
            0.37 3];
        
        root_mult_arr_gx = [
            0.56    20
            0.75    3
            0.99    4
            0.37    3
            0.12    3
            0.20    3
            ];
        
    case '7'
        
        root_mult_arr_fx = [0.1    20
            0.5    2];
        
        root_mult_arr_gx = [0.1    20
            0.9    1];
        
    case '8'
        
        root_mult_arr_fx = [0.1    2
            0.3     2
            0.5    2];
        
        root_mult_arr_gx = [0.1    2];
        
        % From Winkler Paper - Methods for the computation of the degree of an
        % approximate greatest common divisor of two inexact bernstein basis
        % polynomials.
    case '9'
        root_mult_arr_fx = [
            0.10    3
            0.56    4
            0.75    3
            0.82    3
            1.37    3
            -0.27   3
            1.46    2
            ];
        
        root_mult_arr_gx = [
            0.10    2
            0.56    4
            0.75    3
            0.99    4
            1.37    3
            2.12    3
            1.20    3
            ];
        
    case '10'
        root_mult_arr_fx = [
            0.23   4
            0.43   3
            0.57   3
            0.92   3
            1.70   3
            ];
        root_mult_arr_gx = [
            0.23   4
            0.30   2
            0.77   5
            0.92   2
            1.20   5
            ];
        
    case '11'
        root_mult_arr_fx = [0.1  5
            0.56 4
            0.75 3
            0.82 3
            1.37 3];
        root_mult_arr_gx = [0.1  5
            0.56 4
            0.75 3
            0.99 4
            1.37 3
            2.12 3
            1.20 3];
        
        % Example 6.2
    case '12'
        root_mult_arr_fx = [
            0.14    3
            0.56    3
            0.89    4
            1.45    4
            2.37    3
            -3.61   5
            ];
        root_mult_arr_gx = [
            0.14    4
            0.99    1
            2.37    3
            -0.76   2
            -1.24   2
            -3.61   7
            ];
        
    case '13'
        root_mult_arr_fx = [
            0.14    3
            0.56    3
            0.89    4
            0.37    3
            ];
        
        root_mult_arr_gx = [
            0.14    3
            0.37    3
            0.76   2
            0.24   2
            ];
        
        
        % Example 6.2
    case '14'
        root_mult_arr_fx = [
            0.10    3
            0.56    5
            1.40    4
            1.79    3
            2.69    2
            -2.68   3
            ];
        root_mult_arr_gx = [
            0.10    4
            0.56    4
            1.79    3
            2.68    2
            2.69    3
            -1.40   2
            ];
    case '15'
        root_mult_arr_fx = [
            0.10    1;
            0.50    1;
            ];
        
        root_mult_arr_gx = [
            ];
        
        
    case 'Custom'
        intvl_low = -1;
        intvl_high = 1;
        
        prompt = 'Enter the degree of Polynomial f(x) :';
        m = input(prompt);
        
        prompt = 'Enter the degree of Polynomial g(x) :';
        n = input(prompt);
        
        prompt = 'Enter the degree of Polynomial d(x) :';
        t = input(prompt);
        
        
        
        [root_mult_arr_fx,root_mult_arr_gx] = BuildRandomPolynomials(m,n,t,intvl_low, intvl_high);
        
        
end

roots_dx = GetDivisor(root_mult_arr_fx,root_mult_arr_gx);
roots_ux = GetQuotient(root_mult_arr_fx,roots_dx);
roots_vx = GetQuotient(root_mult_arr_gx,roots_dx);

m = sum(root_mult_arr_fx(:,2));
n = sum(root_mult_arr_gx(:,2));
d = sum(roots_dx(:,2));
rankloss = d;

PrintFactorization(root_mult_arr_fx,'f');
PrintFactorization(root_mult_arr_gx,'g');
PrintFactorization(roots_dx,'d');

end




