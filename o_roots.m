function [] = o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_roots(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
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
%
% % Example
% >> o_roots('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', 'y', 'Standard STLN')
% >> o_roots('Custom:m=5 low=-1 high=2', 1e-10, 1e-12, 'Geometric Mean Matlab Method', 'y', 'Standard STLN')

% Add Subfolders
addpath('Root Finding',...
    'Deconvolution',...
    'Examples',...
    'Formatting',...
    'GCD Finding',...
    'LowRankApproximation');

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
SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method)

% %
% %
% Get polynomial f(x)
fx = Examples_Roots(ex_num);

% Add Noise to f(x)
fx = Noise(fx,emin);

% %
% % MY METHOD
% %
% %

% try
    mymethod_tic = tic;
    % Get roots by my method and compare the computed f(x) with the exact f(x)
    [root_multiplicity_array_mymethod] = o_roots_mymethod(fx);
    rel_err.MyMethod = GetRelativeError(root_multiplicity_array_mymethod,fx,'My Method');
    LineBreakLarge()
    time.MyMethod = toc(mymethod_tic);
% catch err
%     fprintf([err.message '\n'])
%     fprintf('Error my method \n')
%     rel_err.MyMethod = 9999999;
%     time.MyMethod = 9999999;
% end


% %
% % MUSSER METHOD
% %
% %

try
    MusserMethod_tic = tic;
    [root_multiplicity_array_MusserMethod] = o_roots_Musser(fx);
    rel_err.MusserMethod = GetRelativeError(root_multiplicity_array_MusserMethod,fx,'Musser Method');
    LineBreakLarge()
    time.MusserMethod = toc(MusserMethod_tic);
catch err
    
    fprintf([mfilename ' : ' sprintf('Error in Musser Method \n')])
    fprintf(err.message);
    rel_err.MusserMethod = 999999;
    time.MusserMethod = 999999;
    
end

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

fullFileName = 'Results_o_roots.txt';


if exist('Results_o_roots.txt', 'file')
    
    fileID = fopen('Results_o_roots.txt','a');
    
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
        SETTINGS.DECONVOLUTION_METHOD,...
        time.MyMethod,...
        time.MusserMethod,...
        SETTINGS.ROOTS_UX...
        );
    fclose(fileID);
    
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end
