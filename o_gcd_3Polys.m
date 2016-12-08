function [] = o_gcd_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method)
% o_gcd_3Polys(ex_num,el,mean_method, bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Given THREE polynomials f(x) and g(x) calculate the GCD d(x).
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
% >> o_gcd_3Polys('1',1e-12,1e-10,'Geometric Mean Matlab Method', 'y','Standard STLN','Standard APF Nonlinear')
%
% >> o_gcd_3Polys(ex_num,1e-12,1e-10,'Geometric Mean Matlab Method', 'y','Standard STLN', 'Standard APF Nonlinear')
% >> ex_num = 'Custom:m=10n=5t=2.low=-1high=2'

% Initialise global variables
global SETTINGS

% Set Problem Type : Either 'GCD' or 'Roots'
problem_type = 'GCD 3 Polynomials';

% % Add subfolders
restoredefaultpath

addpath('Build Matrices',...
    'Formatting',...
    'Get Cofactor Coefficients',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...
    'GCD Finding',...
    'Plotting',...
    'Preprocessing');

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Low Rank Approximation'));

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)

% Get coefficients of f(x,y) g(x,y) from example file
[fx_exact, gx_exact, hx_exact, dx_exact, ux_exact, vx_exact, wx_exact] = Examples_GCD_FromCoefficients_3Polys(ex_num);



% Add noise to the coefficients of polynomials f(x) and g(x) at a
% predefined signal-to-noise ratio.
[fx_noisy,~] = AddVariableNoiseToPoly(fx_exact,emin,emax);
[gx_noisy,~] = AddVariableNoiseToPoly(gx_exact,emin,emax);
[hx_noisy,~] = AddVariableNoiseToPoly(hx_exact,emin,emax);

fx = fx_noisy;
gx = gx_noisy;
hx = hx_noisy;

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
upper_bound = min([GetDegree(fx),GetDegree(gx),GetDegree(hx)]);
lower_bound = 1;
deg_limits = [lower_bound,upper_bound];

% Compute degree of gcd by my method
[fx_calc, gx_calc, hx_calc, dx_calc, ux_calc, vx_calc, wx_calc, ~, ~] = o_gcd_mymethod_3Polys(fx,gx,hx,deg_limits);

% Get distance of the computed d(x) from the exact d(x)

my_error.dx = GetDistance(dx_exact,dx_calc);
my_error.ux = GetDistance(ux_exact,ux_calc);
my_error.vx = GetDistance(vx_exact,vx_calc);
my_error.wx = GetDistance(wx_exact,wx_calc);

my_error.MyMethod = my_error.dx;
display(my_error); 


PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(dx_exact),GetDegree(dx_calc),my_error)

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
        my_error('PLOT_GRAPHS is either y or n');
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
try
dist = norm(f_exact - f_computed) ./ norm(f_exact);
catch
    dist = 1000;
end
end


function [] = PrintToFile(m,n,t_exact,t_comp,error)
%
%
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Computed degree of the GCD d(x)
%
% error : array of errors e
%   error.dx
%   error.ux
%   error.vx

global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd%s.txt',datetime('today'));

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end

    function WriteNewLine()
        
        % 19 fields
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...,...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(t_exact),...    
            int2str(t_comp),...
            num2str(error.ux),...
            num2str(error.vx),...
            num2str(error.dx),...
            SETTINGS.MEAN_METHOD,...
            SETTINGS.BOOL_ALPHA_THETA,...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.BOOL_LOG,...
            SETTINGS.GCD_COEFFICIENT_METHOD...
            );
            
    end

    function WriteHeader()
        fprintf(fileID,'DATE,EX_NUM,m,n,t_exact,t_comp,ERROR_UX,ERROR_VX,ERROR_DX,MEAN_METHOD,BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,BOOL_LOG,GCD_METHOD \n');
    end







end

