function [fx,gx,dx,ux,vx] = Examples_GCD(ex_num)
% Given an example number, return the coefficients of two polynomials

% Output :
%
% fx : Column vector where row i contains coefficient a_{i} x^{i}

root_mult_array_fx = [];
root_mult_array_gx = [];
root_mult_array_dx = [];
root_mult_array_ux = [];
root_mult_array_vx = [];

switch ex_num
    case '1'
        root_mult_array_fx = ...
            [
            1 2;
            1.2 1;
            ];
        root_mult_array_gx = ...
            [
            1 2;
            2   1];
        
    case '2'
        root_mult_array_fx = ...
            [
            0.5 4;
            1.2 3;
            0.1564  2;
            0.78615871145   3
            ];
        root_mult_array_gx = ...
            [
            0.5 4;
            2   2
            0.1564  2;
            0.78615871145   3
            ];
        
        
    case 'Sederberg'
        root_mult_array_fx = ...
            [
            6   1;
            -3  1;
            -7  1;
            ];
        root_mult_array_gx = ...
            [
            15  1;
            -5  1;
            -10 1;
            ];
    otherwise
        error('not a valid example number')
end

root_mult_array_dx = GetDivisor(root_mult_array_fx,root_mult_array_gx);
root_mult_array_ux = GetQuotient(root_mult_array_fx,root_mult_array_dx);
root_mult_array_vx = GetQuotient(root_mult_array_gx,root_mult_array_dx);

% Given the roots and multiplicities of f(x), get the coefficients
fx = GetCoefficients(root_mult_array_fx);
gx = GetCoefficients(root_mult_array_gx);
dx = GetCoefficients(root_mult_array_dx);
ux = GetCoefficients(root_mult_array_ux);
vx = GetCoefficients(root_mult_array_vx);

end

