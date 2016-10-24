function [fx, gx, dx,ux,vx] = Examples_GCD(ex_num)
% Get a set of polynomials f(x),g(x),d(x),u(x),v(x) given an example number
%
% Input :
%
% ex_num : (String) - Note that the example number may be non-numeric.
% 'Custom' allows user to define the degree of f(x) g(x) and d(x).
%
%
% Output :
%
% where row i contains coefficient a_{i} x^{i}
%
% fx : Column vector of coefficients of polynomial f(x)
%
% gx : Column vector of coefficients of polynomial g(x)
%
% dx : Column vector of coefficients of polynomial d(x) where d(x) is the
% GCD of f(x) and g(x).
%
% ux : Column vector of coefficients of polynomial u(x) where u(x) is given
% by f(x) divided by d(x).
%
% vx : column vector of coefficients of polynomail v(x) where v(x) is given
% by g(x) divided by d(x)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        % Get inputs f and g and
        [...
            root_mult_array_fx,...
            root_mult_array_gx,...
            root_mult_array_dx,...
            root_mult_array_ux,...
            root_mult_array_vx] = Examples_GCD_FromRoots(ex_num);
        
               
        % Given the roots and multiplicities of f(x), get the coefficients
        fx = GetCoefficientsFromRoots(root_mult_array_fx);
        gx = GetCoefficientsFromRoots(root_mult_array_gx);
        dx = GetCoefficientsFromRoots(root_mult_array_dx);
        ux = GetCoefficientsFromRoots(root_mult_array_ux);
        vx = GetCoefficientsFromRoots(root_mult_array_vx);

        
        
    case 'From Coefficients'
        
        [fx,gx,dx] = Examples_GCD_FromCoefficients(ex_num);
        ux = 1;
        vx = 1;
        
        
    otherwise
        error('err')
        
end

end



