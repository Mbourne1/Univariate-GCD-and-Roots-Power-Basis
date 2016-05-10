function [] = o_gcd(ex_num,el,mean_method, bool_alpha_theta,low_rank_approx_method)
% o_gcd(ex_num,el,mean_method, bool_alpha_theta,low_rank_approx_method)
%
% Given two polynomials f(x) and g(x) calculate the GCD d(x).
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% el : Noise lower threshold
%
% mean_method :
%       'None'
%       'Geometric Mean Matlab Method'
%       'Geometric Mean My Method'
%
% bool_alpha_theta : 'y' or 'n' (Include/ Exclude Preprocessing)
%
% low_rank_approx_method:
%       'Standard SNTLN'
%       'Standard STLN'
%       'None'
%
%
% % Example
%
% >> o_gcd('1',1e-12,'Geometric Mean Matlab Method', 'y','Standard STLN')

global SETTINGS

SETTINGS.NOISE = el;
SETTINGS.EX_NUM= ex_num;

SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)

% Get the example
[fx_exact,gx_exact,dx_exact,ux_exact,vx_exact] = Examples_GCD(ex_num);


% Add noise to the coefficients of polynomials f(x) and g(x) at a
% predefined signal-to-noise ratio.
[fx,~] = Noise(fx_exact,el);
[gx,~] = Noise(gx_exact,el);

% %
% %
% %
% Get the GCD by Zeng method
%[d_Zeng,v_Zeng,w_Zeng,res,cond] = o_gcd_zeng(fx,gx);

% %
% %
% %
% Get the GCD d(x) of f(x) and g(x) by my method

% Get upper and lower bound of degree of GCD.
upper_bound = min(GetDegree(fx),GetDegree(gx));
lower_bound = 1;
deg_limits = [lower_bound,upper_bound];

% Compute degree of gcd by my method
[~,~,dx_calc,ux_calc,vx_calc,~,~] = o_gcd_mymethod(fx,gx,deg_limits);

% %
% %
% %
% Print the coefficients of d(x), u(x) and v(x)
PrintCoefficients(ux_exact,ux_calc,'u(x)')
GetDistance(ux_exact,ux_calc,'u(x)');

%
PrintCoefficients(vx_exact,vx_calc,'v(x)')
GetDistance(vx_exact,vx_calc,'v(x)');

%
PrintCoefficients(dx_exact,dx_calc,'d(x)')
GetDistance(dx_exact,dx_calc,'d(x)');

%
%GetDistance(dx_exact,d_Zeng,'d(x)');

% Get distance of the computed d(x) from the exact d(x)
error_dx = GetDistance(dx_exact,dx_calc,'d(x)');

PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(dx_calc),error_dx)


%% Plot the three curves f(x), g(x) and d(x)

switch SETTINGS.PLOT_GRAPHS
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

function [dist] = GetDistance(f_exact,f_computed)
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


function [] = PrintToFile(m,n,t,error_dx)

global SETTINGS

fullFileName = 'o_gcd_results.txt';


if exist('o_gcd_results.txt', 'file')
    fileID = fopen('o_gcd_results.txt','a');
    fprintf(fileID,'%s \t %5d \t %5d \t %5d \t %s \t %s \t %s \t %s \t %s \n',...
        SETTINGS.EX_NUM,...
        m,n,t,error_dx,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.NOISE,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
    fclose(fileID);

else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end

