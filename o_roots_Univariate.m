function [] = o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
% O_ROOTS_UNIVARIATE(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
%
% Calculate the roots of a polynomial f(x) by a series of methods.
%
% % Inputs.
%
% ex_num : Example Number
%          Note a cusom example can be specified using the reg exp
%          'Custom:m=(\d+).low=(-?\d+).high=(-?\d+)'
%   
% emin : Lower Noise Level
%
% emax : Upper Noise Level
%
% mean_method : Method of mean used to divide coefficients of Sylvester
%               Matrix as part of preprocessing.
%   'None' : No mean method is used
%   'Geometric Mean Matlab Method' : Divide by Geometric Mean
%
% bool_alpha_theta :
%   'y' : Include Preprocessing
%   'n' : Exclude Preprocessing
%
% low_rank_approx_method :
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
%   'Root Specific SNTLN'
%
% apf_method
%   'None'
%   'Standard APF Nonlinear' - Not Developed
%   'Standard APF Linear' - Not Developed
%
% % Example
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'None', false, 'None','None')
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None','None')
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN','Standard APF Nonlinear')
% >> O_ROOTS_UNIVARIATE('Custom:m=5 low=-1 high=2', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN','Standard APF Nonlinear')

% Add Subfolders
restoredefaultpath
addpath(...
    'Build Matrices',...
    'Deconvolution',...
    'Formatting',...
    'GCD Finding',...
    'Get Cofactor Coefficients',...
    'Get GCD Coefficients',...
    'Low Rank Approximation',...
    'Plotting',...
    'Preprocessing'...
    );

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Get GCD Degree'));
addpath(genpath('Root Finding'));
addpath(genpath('Low Rank Approximation'));

% Initialise global variables
global SETTINGS

if ~isempty(SETTINGS)
else
    % If global variables are not yet set
    SETTINGS.DECONVOLUTION_METHOD = 'Batch';
    SETTINGS.ROOTS_UX = 'From ux';
end

% Set the problem type to be of type 'ROOTS'
problem_type = 'Roots';

% Set global input by user.
SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)

% %
% %
% Get polynomial f(x)
fx = Examples_Roots(ex_num);

% Add noise to coefficients of f(x)
fx = AddNoiseToPoly(fx,emin);

% %
% % MY METHOD
% %
% %

 %try
    mymethod_tic = tic;
%    Get roots by my method and compare the computed f(x) with the exact f(x)
    [root_multiplicity_array_mymethod] = o_roots_mymethod(fx);
    rel_err.MyMethod = GetRelativeError(root_multiplicity_array_mymethod,fx,'My Method');
    LineBreakLarge()
    time.MyMethod = toc(mymethod_tic);
%  catch err
%      fprintf([err.message '\n'])
%      fprintf('Error my method \n')
%      rel_err.MyMethod = 9999999;
%      time.MyMethod = 9999999;
%  end


% %
% % MUSSER METHOD
% %
% %

%try
    MusserMethod_tic = tic;
    [root_multiplicity_array_MusserMethod] = o_roots_Musser(fx);
    rel_err.MusserMethod = GetRelativeError(root_multiplicity_array_MusserMethod,fx,'Musser Method');
    LineBreakLarge()
    time.MusserMethod = toc(MusserMethod_tic);
%catch err
   
%    fprintf([mfilename ' : ' sprintf('Error in Musser Method \n')])
%    fprintf([err.message '\n']);
%    rel_err.MusserMethod = 999999;
%    time.MusserMethod = 999999;
    
%end

% %
% % YUN METHOD
% %
% %

% YunMethod_tic = tic;
% [root_multiplicity_array_YunMethod] = o_roots_Yun(fx);
% rel_err.YunMethod = GetRelativeError(root_multiplicity_array_YunMethod,fx_exact,'Musser Method');
% LineBreakLarge()
% YunMethod_duration = toc(YunMethod_tic);


LineBreakLarge();
fprintf('Duration - My Method : %2.4f\n', time.MyMethod);
fprintf('Duration - Musser Method : %2.4f\n', time.MusserMethod);
LineBreakLarge();


% %
% % MATLAB
% Get roots by matlab method
root_mult_array_MatlabMethod = o_roots_matlab(fx);
rel_err.MatlabMethod = GetRelativeError(root_mult_array_MatlabMethod,fx,'Matlab Method');
LineBreakLarge()


% %
% % MULTROOT
% Get roots by zheng method
computed_root_mult_array_multroot = o_roots_multroot(flipud(fx));
rel_err.MultrootMethod = GetRelativeError(computed_root_mult_array_multroot,fx,'MultRoot Method');
LineBreakLarge()


PrintToFile(rel_err,time)

%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
if( SETTINGS.PLOT_GRAPHS)
    
        figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
        figure('name',figure_name)
        hold on;
        scatter( real(root_multiplicity_array_mymethod(:,1)), imag(root_multiplicity_array_mymethod(:,1)),'yellow','*','DisplayName','My Method');
        scatter( real(root_mult_array_MatlabMethod(:,1)), imag(root_mult_array_MatlabMethod(:,1)),'red','DisplayName','Matlab Roots');
        scatter( real(computed_root_mult_array_multroot(:,1)), imag(computed_root_mult_array_multroot(:,1)),'green','s','filled','DisplayName','MultRoots');
        
       
        xlabel('Real');
        ylabel('Imaginary');
        legend(gca,'show')
        ylim()
        str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
        title(str);
        hold off
    
end


end

function [rel_err] = GetRelativeError(computed_root_mult_array,fx_exact,method)
% Get distance between the computed polynomial, and the exact polynomial,
% given the computed roots.
%
% Inputs.
%
% computed_root_mult_array
%
% method : (String)

try
    fx_computed = GetCoefficients(computed_root_mult_array);
    
    fprintf('Distance between f_exact and f_comp by %s : \n',method)
    
    % Get the relative error.
    rel_err = norm(fx_exact - fx_computed) ./ norm(fx_exact) ;
    
    display(rel_err);
catch
    rel_err = 100000000000000;
end
end


function [] = PrintToFile(rel_err,time)

global SETTINGS
addpath 'Results'

fullFileName = 'Results/Results_o_roots.txt';


if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    format_str = ...
        '%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n';
    fprintf(fileID,format_str,...
        datetime('now'),...
        SETTINGS.PROBLEM_TYPE,...
        SETTINGS.EX_NUM,...
        rel_err.MyMethod,...
        rel_err.MusserMethod,...
        rel_err.MatlabMethod,...
        rel_err.MultrootMethod,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        SETTINGS.DECONVOLUTION_METHOD_FX_HX,...
        time.MyMethod,...
        time.MusserMethod,...
        SETTINGS.ROOTS_HX...
        );
    fclose(fileID);
    
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end
