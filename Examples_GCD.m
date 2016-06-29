function [noisy_fx, noisy_gx, fx, gx, dx, ux, vx] = Examples_GCD(ex_num)
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

EXAMPLE_TYPE = 'FromRoots';
global SETTINGS

switch EXAMPLE_TYPE
    case 'FromRoots'
        % Get inputs f and g and
        [...
            root_mult_array_fx,...
            root_mult_array_gx,...
            root_mult_array_dx,...
            root_mult_array_ux,...
            root_mult_array_vx] = Examples_GCD_Roots(ex_num);
        
        noise = SETTINGS.EMIN;
        
        % Given the roots and multiplicities of f(x), get the coefficients
        fx = GetCoefficients(root_mult_array_fx);
        gx = GetCoefficients(root_mult_array_gx);
        dx = GetCoefficients(root_mult_array_dx);
        ux = GetCoefficients(root_mult_array_ux);
        vx = GetCoefficients(root_mult_array_vx);
        
        % Add noise to the roots
        switch SETTINGS.NOISE_TYPE
            
            case 'Roots'
                noisy_root_mult_array_fx = AddNoiseToRoots(root_mult_array_fx,noise);
                noisy_root_mult_array_gx = AddNoiseToRoots(root_mult_array_gx,noise);
                
                noisy_fx = GetCoefficients(noisy_root_mult_array_fx);
                noisy_gx = GetCoefficients(noisy_root_mult_array_gx);
                
            case 'Coefficients'
                
                emin = SETTINGS.EMIN;
                emax = SETTINGS.EMAX;
                
                % Add noise to the coefficients of polynomials f(x) and g(x) at a
                % predefined signal-to-noise ratio.
                [noisy_fx,~] = Noise(fx,emin,emax);
                [noisy_gx,~] = Noise(gx,emin,emax);
                
            otherwise
                error('Noise type either Roots or coefficients');
        end
        
        
        fprintf([mfilename fprintf('dist f and f_noisy')]);
        display(norm(fx-noisy_fx) ./ norm(fx))
        display(norm(gx-noisy_gx) ./ norm(gx))
        
        
    case 'FromCoefficients'
        
        [fx,gx,~] = GCD_Examples_Coefficients(ex_num);
        
        
end

end


function [root_mult_array_fx,root_mult_array_gx,root_mult_array_dx,root_mult_array_ux,root_mult_array_vx] = Examples_GCD_Roots(ex_num)
% Given an example number, return the coefficients of two polynomials
%
% Inputs:
%
% Example Number : (int)
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

pattern = 'Custom:m=(\d+).n=(\d+).t=(\d+).low=(-?\d+).high=(-?\d+)';

if ~isempty(regexp(ex_num,pattern,'start'))
    
    str = ex_num;
    
    expression_m = regexp(str,'m=(\d+)','tokens');
    m_str = expression_m{1};
    
    expression_n = regexp(str,'n=(\d+)','tokens');
    n_str = expression_n{1};
    
    expression_t = regexp(str,'t=(\d+)','tokens');
    t_str = expression_t{1};
    
    expression_low = regexp(str, 'low=(-?\d+)', 'tokens');
    low_str = expression_low{1};
    
    expression_high = regexp(str, 'high=(-?\d+)','tokens');
    high_str = expression_high{1};
    
    
    
    m = str2double(m_str);
    n = str2double(n_str);
    t = str2double(t_str);
    intvl_low = str2double(low_str);
    intvl_high = str2double(high_str);
    
    [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,n,t,intvl_low, intvl_high);
    
    
    
    
else
    
    switch ex_num
        
        case '0'
            root_mult_array_fx = ...
            [
                1.456   1
                0.567   2
                0.927   3
            ];
        
        
            root_mult_array_gx = ...
            [
                1.456   1
                0.567   2
                0.427   3
            ];
        
        
        
                       
            
        case '1' % From Bini
            root_mult_array_fx = ...
                [
                -1      1;
                -2     1;
                2      1;
                1      1;
                ];
            
            root_mult_array_gx = ...
                [
                2      1;
                -0.5   1;
                0.5    1;
                ];
            
        case '2'
            root_mult_array_fx = ...
                [
                -1      1;
                -2     1;
                2      1;
                1      1;
                3      1;
                -3     1;
                ];
            
            root_mult_array_gx = ...
                [
                2      1;
                -0.5   1;
                0.5    1;
                3      1;
                -3     1;
                ];
            
        case '3'
            root_mult_array_fx = ...
                [
                0.10    10
                0.56    4
                0.40    4
                0.79    3
                0.69    2
                ];
            
            root_mult_array_gx = ...
                [
                0.10    10
                0.56    4
                0.69    2
                ];
            
            
            
        case '4'
            root_mult_array_fx = [
                0.2 2
                0.5 1
                0.7 1
                0.9 1
                1.1 1
                ];
            root_mult_array_gx = ...
                [
                0.5 1
                0.1 1
                0.9 1
                0.3 1
                ];
            
            
        case '5'
            root_mult_array_fx = ...
                [
                0.2 1
                0.4 1
                0.6 1
                0.8 1
                ];
            
            root_mult_array_gx = ...
                [
                0.9 1
                0.2 1
                0.3 1
                ];
            
        case '6'
            root_mult_array_fx = ...
                [
                0.56 20
                0.75 3
                0.82 3
                0.37 3
                ];
            
            root_mult_array_gx = ...
                [
                0.56    20
                0.75    3
                0.99    4
                0.37    3
                0.12    3
                0.20    3
                ];
            
        case '7'
            
            root_mult_array_fx = ...
                [
                0.1    20
                0.5    2
                ];
            
            root_mult_array_gx = ...
                [
                0.1    20
                0.9    1
                ];
            
        case '8'
            
            root_mult_array_fx = ...
                [
                0.1    2
                0.3     2
                0.5    2
                ];
            
            root_mult_array_gx = [0.1    2];
            
            % From Winkler Paper - Methods for the computation of the degree of an
            % approximate greatest common divisor of two inexact bernstein basis
            % polynomials.
        case '9'
            root_mult_array_fx = ...
                [
                0.10    3
                0.56    4
                0.75    3
                0.82    3
                1.37    3
                -0.27   3
                1.46    2
                ];
            
            root_mult_array_gx = ...
                [
                0.10    2
                0.56    4
                0.75    3
                0.99    4
                1.37    3
                2.12    3
                1.20    3
                ];
            
        case '10'
            root_mult_array_fx = ...
                [
                0.23   4
                0.43   3
                0.57   3
                0.92   3
                1.70   3
                ];
            root_mult_array_gx = ...
                [
                0.23   4
                0.30   2
                0.77   5
                0.92   2
                1.20   5
                ];
            
        case '11'
            root_mult_array_fx = ...
                [
                0.1  5
                0.56 4
                0.75 3
                0.82 3
                1.37 3
                ];
            root_mult_array_gx = ...
                [
                0.1  5
                0.56 4
                0.75 3
                0.99 4
                1.37 3
                2.12 3
                1.20 3
                ];
            
            % Example 6.2
        case '12'
            root_mult_array_fx = ...
                [
                0.14    3
                0.56    3
                0.89    4
                1.45    4
                2.37    3
                -3.61   5
                ];
            root_mult_array_gx = ...
                [
                0.14    4
                0.99    1
                2.37    3
                -0.76   2
                -1.24   2
                -3.61   7
                ];
            
        case '13'
            root_mult_array_fx = ...
                [
                0.14    3
                0.56    3
                0.89    4
                0.37    3
                ];
            
            root_mult_array_gx = [
                0.14    3
                0.37    3
                0.76   2
                0.24   2
                ];
            
            
            % Example 6.2
        case '14'
            root_mult_array_fx = [
                0.10    3
                0.56    5
                1.40    4
                1.79    3
                2.69    2
                -2.68   3
                ];
            root_mult_array_gx = [
                0.10    4
                0.56    4
                1.79    3
                2.68    2
                2.69    3
                -1.40   2
                ];
        case '15'
            root_mult_array_fx = [
                0.10    1;
                0.50    1;
                ];
            
            root_mult_array_gx = [
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
            
        case 'Custom'
            intvl_low = -5;
            intvl_high = 5;
            
            prompt = 'Enter the degree of Polynomial f(x) :';
            m = input(prompt);
            
            prompt = 'Enter the degree of Polynomial g(x) :';
            n = input(prompt);
            
            prompt = 'Enter the degree of Polynomial d(x) :';
            t = input(prompt);
            
            
            
            [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,n,t,intvl_low, intvl_high);
        otherwise
            error(sprintf([mfilename ' : Not a valid example number' ex_num]))
    end
    
    %PrintFactorization(root_mult_array_fx,'f');
    %PrintFactorization(root_mult_array_gx,'g');
    
end

root_mult_array_dx = GetDivisor(root_mult_array_fx,root_mult_array_gx);
root_mult_array_ux = GetQuotient(root_mult_array_fx,root_mult_array_dx);
root_mult_array_vx = GetQuotient(root_mult_array_gx,root_mult_array_dx);

display(root_mult_array_fx)
display(root_mult_array_gx)
display(root_mult_array_dx)

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

function [fx,gx,dx_exact] = GCD_Examples_Coefficients(ex_num)
switch ex_num
    case '1'
        
        fx = [1; 1; -2; -2; 1; 1];
        gx = [-2; 1; 4; -2; -2; 1];
        dx_exact = [1; 0; -2; 0; 1];
    case '2'
        
        fx = [-1; -1; 1; 1];
        gx = [2; -1; -2; 1];
        dx_exact = [-1; 0; 1];
        
    case '3'
        fx = [-5; 1; -5; 1];
        gx = [-2; 1; -2; 1];
        dx_exact = [1; 0; 1];
        
    case '4'
        fx = [6.8; -17.6; 16.5; -6.7; 1];
        gx = [6;   -16  ; 15.5; -6.5; 1];
        
end

end

