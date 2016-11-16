function [] = o_gcd(ex_num,emin,emax,mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
% o_gcd(ex_num,el,mean_method, bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Given two polynomials f(x) and g(x) calculate the GCD d(x).
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% el : Noise lower level
%
% em : Noise upper level
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
% apf_method : 
%       'Standard Nonlinear'
%       'Standard Linear'
%       'None'
%
% % Example
%
% >> o_gcd('1',1e-12,1e-10,'Geometric Mean Matlab Method', 'y','Standard STLN','Standard APF Nonlinear')
%
% >> o_gcd(ex_num,1e-12,1e-10,'Geometric Mean Matlab Method', 'y','Standard STLN', 'Standard APF Nonlinear')
% >> ex_num = 'Custom:m=10n=5t=2.low=-1high=2'

% Initialise global variables
global SETTINGS

% Set Problem Type : Either 'GCD' or 'Roots'
problem_type = 'GCD';

% Set global variables
SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)



% Add subfolders
restoredefaultpath

% Add Relevant Paths
addpath('Build Matrices',...
    'Formatting',...
    'Get GCD Degree',...
    'GCD Finding',...
    'Plotting',...
    'Preprocessing');

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Low Rank Approximation'));


% Get coefficients of f(x,y) g(x,y) from example file
[fx_exact, gx_exact, dx_exact,ux_exact,vx_exact] = Examples_GCD(ex_num);


 
emin = SETTINGS.EMIN;
emax = SETTINGS.EMAX;

% Add noise to the coefficients of polynomials f(x) and g(x) at a
% predefined signal-to-noise ratio.
[fx_noisy,~] = AddNoiseToPoly(fx_exact,emin);
[gx_noisy,~] = AddNoiseToPoly(gx_exact,emin);

fx = fx_noisy;
gx = gx_noisy;

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
[fx_calc,gx_calc,dx_calc,ux_calc,vx_calc,~,~] = o_gcd_mymethod(fx,gx,deg_limits);

% Get distance of the computed d(x) from the exact d(x)
try
    
    error_dx = GetDistance(dx_exact,dx_calc);
    error.MyMethod = error_dx;
    error_ux = GetDistance(ux_exact,ux_calc);
    error_vx = GetDistance(vx_exact,vx_calc);
    
    display(error_ux);
    display(error_vx);
    display(error_dx);
    
catch err
    error.MyMethod = 9999999;
    fprintf([mfilename ' : ' err.message '\n'])
end

PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(dx_calc),error)

% %
% %
% %
% Get the GCD by an alternative method
%[dx] = o_gcd_experiment_method(fx,gx)


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

end


function [] = PrintToFile(m,n,t,error)

global SETTINGS

fullFileName = 'Results/Results_o_gcd.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %5s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n',...
        datetime('now'),...,...
        SETTINGS.PROBLEM_TYPE,...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(t),...
        error.MyMethod,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
    fclose(fileID);

else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end

