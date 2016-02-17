function [] = o_gcd(ex_num,preproc,bool_sntln,el)
% Given two polynomials f(x) and g(x) calculate the gcd d(x).


global BOOL_PREPROC
BOOL_PREPROC = preproc;

global BOOL_SNTLN
BOOL_SNTLN = bool_sntln;

% Initialise the Seed for random noise generation
global SEED
SEED = 1024;

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

global BOOL_NOISE
BOOL_NOISE = 'y';

global MAX_ERROR_SNTLN
global MAX_ITE_SNTLN

MAX_ERROR_SNTLN = 1e-12;
MAX_ITE_SNTLN = 100;

EXAMPLE_TYPE = 'FromRoots';

switch EXAMPLE_TYPE
    case 'FromRoots'
        % Get inputs f and g and
        [roots_fx,roots_gx,roots_dx] = GCD_Examples(ex_num);
        
        fx = get_Coeff(roots_fx);
        gx = get_Coeff(roots_gx);
        dx = get_Coeff(roots_dx);
        
    case 'FromCoefficients'
        switch ex_num
            case '1'
                
                fx = [1; 1; -2; -2; 1; 1];
                gx = [-2; 1; 4; -2; -2; 1];
                dx = [1; 0; -2; 0; 1];
            case '2'
                
                fx = [-1; -1; 1; 1]
                gx = [2; -1; -2; 1]
                dx = [-1; 0; 1]
                
            case '3'
                fx = [-5; 1; -5; 1];
                gx = [-2; 1; -2; 1];
                dx = [1; 0; 1];
        end
        
end
% Get degree of polynomial f
[r,~] = size(fx);
m = r - 1;

% Get degree of polynomial g
[r,~] = size(gx);
n = r - 1;

% Add Noise to the coefficients of f(x) and g(x)

switch BOOL_NOISE
    case 'y'
        [fx,f_noise] = noise(fx,el);
        [gx,g_noise] = noise(gx,el);
        
        switch PLOT_GRAPHS
            case 'y'
                
                t = linspace(-100,100,200);
                f_y_noisy = polyval(fx,t);
                f_y_exact = polyval(fx,t);
                
                figure('name','Plotting f vs noisy f')
                hold on
                plot(t,f_y_noisy,'blue')
                plot(t,f_y_exact,'red')
                hold off
                
            case 'n'
        end
    case 'n'
    otherwise
        error('bool_noise is either y or n')
end

% Get the GCD d(x) of f(x) and g(x)
[fx_calc,gx_calc,dx_calc] = o1(fx,gx);

fprintf('Calculated Coefficients of d(x): \n')
dx_calc = dx_calc./dx_calc(1)
fprintf('Exact Coefficients of d(x): \n')
dx./dx(1)


% plot f(x) and g(x)
t = linspace(-10,10,200);
f_y = polyval(fx,t);
g_y = polyval(gx,t);
d_y = polyval(dx_calc,t);

switch PLOT_GRAPHS
    case 'y'
        figure('name','Curve Plot')
        hold on
        plot(t,f_y,'DisplayName','f(y)')
        plot(t,g_y,'DisplayName','g(y)')
        plot(t,d_y,'DisplayName','d(y)')
        legend(gca,'show');
        hold off
    case 'n'
    otherwise
        error('PLOT_GRAPHS is either y or n')
end




end