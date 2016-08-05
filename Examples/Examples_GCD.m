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

EXAMPLE_TYPE = 'FromCoefficients';

switch EXAMPLE_TYPE
    case 'FromRoots'
        % Get inputs f and g and
        [...
            root_mult_array_fx,...
            root_mult_array_gx,...
            root_mult_array_dx,...
            root_mult_array_ux,...
            root_mult_array_vx] = Examples_GCD_FromRoots(ex_num);
        
               
        % Given the roots and multiplicities of f(x), get the coefficients
        fx = GetCoefficients(root_mult_array_fx);
        gx = GetCoefficients(root_mult_array_gx);
        dx = GetCoefficients(root_mult_array_dx);
        ux = GetCoefficients(root_mult_array_ux);
        vx = GetCoefficients(root_mult_array_vx);

        
        
    case 'FromCoefficients'
        
        [fx,gx,dx] = Examples_GCD_FromCoefficients(ex_num);
        ux = 1;
        vx = 1;
        
        
    otherwise
        error('err')
        
end

end


function root_mult_array_fx_noisy = AddNoiseToRoots(root_mult_array_fx, noise)
% Given a set of roots and their corresponding multiplicities. Add noise to
% the roots

roots = root_mult_array_fx(:,1);
mults = root_mult_array_fx(:,2);

nRoots = size(roots,1);

% Generate a vector of random numbers between 1 and -1
a = -1;
b = 1;

rand_vec = (b-a).*rand(nRoots,1) + a;

noise_vec = (roots .* (rand_vec * noise));

noisy_roots = roots + noise_vec;

root_mult_array_fx_noisy = [noisy_roots mults];


end

