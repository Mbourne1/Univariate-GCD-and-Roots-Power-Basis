
function [roots_fx,roots_gx,roots_dx,rankloss,roots_ux,roots_vx]=GCD_Examples(n)
% Inputs.
% n - Index of example to be used
% Outputs.
% a - Roots and multiplicities of polynomial f.
% b - Roots and multiplicities of polynomial g.
% c - Roots and multiplicities of polynomial d, the GCD of f and g.
% u - Roots and multiplicities of quotient polynomial f/u = d.
% v - Roots and multiplicities of quotient polynomial g/v = d.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global SEED

switch n
    
    case '-2' % From Bini
        roots_fx = [
            -1      1;
            -2     1;
            2      1;
            1      1;
            ];
        
        roots_gx = [
            2      1;
            -0.5   1;
            0.5    1;
            ];
    
    case '-1'
        roots_fx = [
            -1      1;
            -2     1;
            2      1;
            1      1;
            3      1;
            -3     1;
            ];
        
        roots_gx = [
            2      1;
            -0.5   1;
            0.5    1;
            3      1;
            -3     1;
            ];

    case '0'
        roots_fx = [
            0.10    30
            0.56    4
            0.40    4
            0.79    3
            0.69    2
            ];
        
        roots_gx = [
            0.10    40
            0.56    4
            0.69    2
            ];

        
        
    case '1'
        roots_fx = [
            0.2 2
            0.5 1
            0.7 1
            0.9 1
            1.1 1
            ];
        roots_gx = [
            0.5 1
            0.1 1
            0.9 1
            0.3 1
            ];

        
    case '2'
        roots_fx = [
            0.2 1
            0.4 1
            0.6 1
            0.8 1
            ];
        
        roots_gx = [0.2 1
            0.3 1];
               
    case '3'
        roots_fx = [
            0.56 20
            0.75 3
            0.82 3
            0.37 3];
        
        roots_gx = [
            0.56    20
            0.75    3
            0.99    4
            0.37    3
            0.12    3
            0.20    3
            ];
     
    case '4'
        
        roots_fx = [0.1    20
            0.5    2];
        
        roots_gx = [0.1    20
            0.9    1];
        
    case '5'
        
        roots_fx = [0.1    2
            0.3     2
            0.5    2];
        
        roots_gx = [0.1    2];
        
        % From Winkler Paper - Methods for the computation of the degree of an
        % approximate greatest common divisor of two inexact bernstein basis
        % polynomials.
    case '6'
        roots_fx = [
            0.10    3
            0.56    4
            0.75    3
            0.82    3
            1.37    3
            -0.27   3
            1.46    2
            ];
        
        roots_gx = [
            0.10    2
            0.56    4
            0.75    3
            0.99    4
            1.37    3
            2.12    3
            1.20    3
            ];
              
    case '7'
        roots_fx = [
            0.23   4
            0.43   3
            0.57   3
            0.92   3
            1.70   3
            ];
        roots_gx = [
            0.23   4
            0.30   2
            0.77   5
            0.92   2
            1.20   5
            ];
         
    case '8'
        roots_fx = [0.1  5
            0.56 4
            0.75 3
            0.82 3
            1.37 3];
        roots_gx = [0.1  5
            0.56 4
            0.75 3
            0.99 4
            1.37 3
            2.12 3
            1.20 3];

    % Example 6.2
    case '9'
        roots_fx = [
            0.14    3
            0.56    3
            0.89    4
            1.45    4
            2.37    3
            -3.61   5
            ];
        roots_gx = [
            0.14    4
            0.99    1
            2.37    3
            -0.76   2
            -1.24   2
            -3.61   7
            ];

    case '10'
        roots_fx = [
            0.14    3
            0.56    3
            0.89    4
            0.37    3
            ];
        
        roots_gx = [
            0.14    3
            0.37    3
            0.76   2
            0.24   2
            ];
        
        
    % Example 6.2
    case '11'
        roots_fx = [
            0.10    3
            0.56    5
            1.40    4
            1.79    3
            2.69    2
            -2.68   3
            ];
        roots_gx = [
            0.10    4
            0.56    4
            1.79    3
            2.68    2
            2.69    3
            -1.40   2
            ];
    case '12'
        roots_fx = [
            0.10    1;
            0.50    1;
        ];
    
        roots_gx = [
            ];
        
    case '999'
        % Degree of f(x)
        m = 6;
        
        % Degree of g(x)
        n = 5;
        
        % Degree of GCD
        t = 3;
        
        % Upper and lower bound of roots
        interval_low    = -10;
        interval_high   =  10;
        
        % Get roots and multiplicities
        [roots_fx,roots_gx] = ...
            BuildRandomPolynomials(m,n,t,interval_low,interval_high,SEED);
        
        
        
end

roots_dx = GetDivisor(roots_fx,roots_gx);

roots_ux = GetQuotient(roots_fx,roots_dx);

roots_vx = GetQuotient(roots_gx,roots_dx);

m = sum(roots_fx(:,2));

n = sum(roots_gx(:,2));

d = sum(roots_dx(:,2));
rankloss = d;

PrintFactorization(roots_fx,'f');
PrintFactorization(roots_gx,'g');
PrintFactorization(roots_dx,'d');

end




