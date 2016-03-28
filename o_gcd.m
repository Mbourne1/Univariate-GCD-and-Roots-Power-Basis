function [] = o_gcd(ex_num,el,bool_preproc,low_rank_approx_method)
% Given two polynomials f(x) and g(x) calculate the GCD d(x).
%
% Inputs.
%
% ex_num : Example Number
%
% el : Noise lower threshold
%
% bool_preproc : 'y' or 'n' (Include/ Exclude Preprocessing)
%
% low_rank_approx_method: 'Standard SNTLN', 'Standard STLN'
%   

global PLOT_GRAPHS

SetGlobalVariables(bool_preproc,low_rank_approx_method)

EXAMPLE_TYPE = 'FromRoots';

switch EXAMPLE_TYPE
    case 'FromRoots'
        % Get inputs f and g and
        [roots_fx,roots_gx,roots_dx,roots_ux,roots_vx] = GCD_Examples(ex_num);
        
        fx = GetCoefficients(roots_fx);
        gx = GetCoefficients(roots_gx);
        dx_exact = GetCoefficients(roots_dx);
        ux_exact = GetCoefficients(roots_ux);
        vx_exact = GetCoefficients(roots_vx);
        
    case 'FromCoefficients'
        switch ex_num
            case '1'
                
                fx = [1; 1; -2; -2; 1; 1];
                gx = [-2; 1; 4; -2; -2; 1];
                dx_exact = [1; 0; -2; 0; 1];
            case '2'
                
                fx = [-1; -1; 1; 1]
                gx = [2; -1; -2; 1]
                dx_exact = [-1; 0; 1]
                
            case '3'
                fx = [-5; 1; -5; 1];
                gx = [-2; 1; -2; 1];
                dx_exact = [1; 0; 1];
                
            case '4'
                fx = [6.8; -17.6; 16.5; -6.7; 1]
                gx = [6;   -16  ; 15.5; -6.5; 1]
                
        end
        
end

% Get degree of polynomial f
m = GetDegree(fx);

% Get degree of polynomial g
n = GetDegree(gx);

% Add Noise to the coefficients of f(x) and g(x)



% Add noise at a signal to noise ratio.
[fx,~] = Noise(fx,el);
[gx,~] = Noise(gx,el);

% Get the GCD d(x) of f(x) and g(x)
[fx_calc,gx_calc,dx_calc,ux_calc,vx_calc,~,~] = o1(fx,gx);

% Print the coefficients of d(x),u(x) and v(x)


PrintCoefficients(ux_exact,ux_calc,'u(x)')
PrintCoefficients(vx_exact,vx_calc,'v(x)')
PrintCoefficients(dx_exact,dx_calc,'d(x)')

%% Plot the three curves f(x), g(x) and d(x)

switch PLOT_GRAPHS
    case 'y'
        % plot f(x) and g(x)
        t = linspace(-10,10,200);
        f_y = polyval(fx,t);
        g_y = polyval(gx,t);
        d_y = polyval(dx_calc,t);
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

function [] = PrintCoefficients(u_exact,u_computed,name)

fprintf('\nCoefficients of %s \n\n',name);
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')

u_exact = Normalise(u_exact);
u_computed = Normalise(u_computed);

mat = [real(u_exact(:,1))';  real(u_computed(:,1))' ];
fprintf('%f \t \t \t \t %f   \t \t \n', mat);
fprintf('\n');

dist = norm(u_exact - u_computed) ./ norm(u_exact);
fprintf('Distance between %s exact and %s computed : %2.4e \n \n',name,name,dist)

end