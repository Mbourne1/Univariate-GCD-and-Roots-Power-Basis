function [fx,gx,dx,ux,vx] = Examples_GCD(ex_num)
% Given an example number, return the coefficients of two polynomials

% Output :
%
% fx : Column vector where row i contains coefficient a_{i} x^{i}

roots_fx = [];
roots_gx = [];
roots_dx = [];
roots_ux = [];
roots_vx = [];

switch ex_num
    case '1'
        roots_fx = ...
            [
            1 2;
            1.2 1;
            ];
        roots_gx = ...
            [
            1 2;
            2   1];
       
    case '2'
        roots_fx = ...
            [
            0.5 4;
            1.2 3;
            0.1564  2;
            0.78615871145   3
            ];
        roots_gx = ...
            [
            0.5 4;
            2   2
            0.1564  2;
            0.78615871145   3
            ];
       
    
    case 'Sederberg'
        roots_fx = ...
            [
                6   1;
                -3  1;
                -7  1;
            ];
        roots_gx = ...
            [
                15  1;
                -5  1;
                -10 1;
            ];
otherwise
        error('not a valid example number')
end

roots_dx = GetDivisor(roots_fx,roots_gx);
roots_ux = GetQuotient(roots_fx,roots_dx);
roots_vx = GetQuotient(roots_gx,roots_dx);

% Given the roots and multiplicities of f(x), get the coefficients
fx = GetCoefficients(roots_fx);
gx = GetCoefficients(roots_gx);
dx = GetCoefficients(roots_dx);
ux = GetCoefficients(roots_ux);
vx = GetCoefficients(roots_vx);

end

