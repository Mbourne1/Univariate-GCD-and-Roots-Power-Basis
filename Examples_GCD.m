function [fx,gx] = Examples_GCD(ex_num)
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
        roots_dx = ...
            [
            1 2;
            ];
        roots_ux = ...
            [
            1.2 1;
            ];
        roots_vx = ...
            [
            2   1;
            ];
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
        roots_dx = ...
            [
            0.5     4;
            0.1564  2;
            0.78615871145   3
            ];
        roots_ux = ...
            [
            1.2 3;
            ];
        roots_vx = ...
            [
            2   2;
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

% Given the roots and multiplicities of f(x), get the coefficients
fx = get_Coeff(roots_fx);
gx = get_Coeff(roots_gx);
dx = get_Coeff(roots_dx);
ux = get_Coeff(roots_ux);
vx = get_Coeff(roots_vx);

end

