function [] = o_gcd(ex_num,el,mean_method, bool_alpha_theta,low_rank_approx_method)
% o_gcd(ex_num,el,mean_method, bool_alpha_theta,low_rank_approx_method)
%
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

global GCD_OR_ROOTS
GCD_OR_ROOTS = 'GCD';

global PLOT_GRAPHS
global NOISE 
NOISE = el;

SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)

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
        
        [fx,gx,~] = GCD_Examples_Coefficients(ex_num);
        
        
end

% Add noise at a signal to noise ratio.
[fx,~] = Noise(fx,el);
[gx,~] = Noise(gx,el);

% Get GCD by Zeng method
%[u,v,w] = o_gcd_zeng(fx,gx);

% Get the GCD d(x) of f(x) and g(x)
[~,~,dx_calc,ux_calc,vx_calc,~,~] = o1(fx,gx,GetDegree(fx));



% Print the coefficients of d(x),u(x) and v(x)
PrintCoefficients(ux_exact,ux_calc,'u(x)')
GetDistance(ux_exact,ux_calc,'u(x)');

PrintCoefficients(vx_exact,vx_calc,'v(x)')
GetDistance(vx_exact,vx_calc,'v(x)');

PrintCoefficients(dx_exact,dx_calc,'d(x)')
GetDistance(dx_exact,dx_calc,'d(x)');

% Get 
error_dx = GetDistance(dx_exact,dx_calc,'d(x)');

PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(dx_calc),error_dx)


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

function [] = PrintCoefficients(f_exact,f_computed,name)
% Print coefficients of f(x) exact and f(x) computed.
% 
% Inputs.
%
% f_exact : Coefficients of polynomial f(x)
%
% f_computed : Coefficients of polynomial f(x) computed
% 
% name : 

fprintf('\nCoefficients of %s \n\n',name);
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')

f_exact = Normalise(f_exact);
f_computed = Normalise(f_computed);

mat = [real(f_exact(:,1))';  real(f_computed(:,1))' ];
fprintf('%f \t \t \t \t %f   \t \t \n', mat);
fprintf('\n');

distance = GetDistance(f_exact,f_computed,name);
fprintf('Distance between %s exact and %s computed : %2.4e \n \n',name,name,distance)


end

function [dist] = GetDistance(f_exact,f_computed,name)
% GetDistance(u_exact,u_computed,name)
%
% Get the distance between the coefficients of two vectors.
%
% Inputs.
%
% f_exact : 
%
% f_computed : 
%
% name : Name of function used for printing

% Normalise f(x) and f(x) computed.
f_exact = Normalise(f_exact);
f_computed = Normalise(f_computed);

% Get distance
dist = norm(f_exact - f_computed) ./ norm(f_exact);

% Print

end


function []= PrintToFile(m,n,t,error_dx)

global NOISE
global BOOL_PREPROC
global LOW_RANK_APPROXIMATION_METHOD

fullFileName = 'o_gcd_results.txt';


if exist('o_gcd_results.txt', 'file')
    fileID = fopen('o_gcd_results.txt','a');
    fprintf(fileID,'%5d \t %5d \t %5d \t %s \t %s \t %s \t %s \n',...
        m,n,t,error_dx,BOOL_PREPROC, NOISE,LOW_RANK_APPROXIMATION_METHOD);
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end




end

