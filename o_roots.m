function [] = o_roots(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_roots(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
%
% Calculate the roots of a polynomial f(x) by a series of methods.
%
% % Inputs.
%
% ex_num : Example Number
%
% el : Lower Noise Level
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
% >> o_roots('1',1e-10,'Geometric Mean Matlab Method', 'y','Standard STLN')
% 


% %
% Set global variables.
SetGlobalVariables(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)

% %
% Get polynomial f(x)
fx_exact = Examples_Roots(ex_num);

% Add Noise to f(x)
fx = Noise(fx_exact,el);

% %
% %
% %
% MY METHOD
mymethod_tic = tic;
% Get roots by my method and compare the computed f(x) with the exact f(x)
[root_multiplicity_array_mymethod] = o_roots_mymethod(fx);
rel_err.MyMethod = GetRelativeError(root_multiplicity_array_mymethod,fx_exact,'My Method');
LineBreakLarge()
mymethod_duration = toc(mymethod_tic);

% %
% %
% %
% MUSSER METHOD
MusserMethod_tic = tic;
[root_multiplicity_array_MusserMethod] = o_roots_Musser(fx);
rel_err.MusserMethod = GetRelativeError(root_multiplicity_array_MusserMethod,fx_exact,'Musser Method');
LineBreakLarge()
MusserMethod_duration = toc(MusserMethod_tic);

LineBreakLarge();
fprintf('Duration - My Method : %2.4f\n', mymethod_duration);
fprintf('Duration - Musser Method : %2.4f\n', MusserMethod_duration);
LineBreakLarge();
% %
% % MATLAB
% Get roots by matlab method

root_mult_array_MatlabMethod = o_roots_matlab(fx);
rel_err.MatlabMethod = GetRelativeError(root_mult_array_MatlabMethod,fx_exact,'Matlab Method');
LineBreakLarge()
% %
% % MULTROOT
% Get roots by zheng method
computed_root_mult_array_multroot = o_roots_multroot(flipud(fx));
rel_err.MultrootMethod = GetRelativeError(computed_root_mult_array_multroot,fx_exact,'MultRoot Method');
LineBreakLarge()


PrintToFile(rel_err)

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


fx_computed = GetCoefficients(computed_root_mult_array);

fprintf('Distance between f_exact and f_comp by %s : \n',method)

% Get the relative error.
rel_err = norm(fx_exact - fx_computed) ./ norm(fx_exact) ;

display(rel_err);
end


function [] = PrintToFile(rel_err)

global SETTINGS

fullFileName = 'o_roots_results.txt';


if exist('o_roots_results.txt', 'file')

    fileID = fopen('o_roots_results.txt','a');

    format_str = ...
        '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n';
    fprintf(fileID,format_str,...
        SETTINGS.EX_NUM,...
        rel_err.MyMethod,...
        rel_err.MusserMethod,...
        rel_err.MatlabMethod,...
        rel_err.MultrootMethod,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.NOISE,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        SETTINGS.DECONVOLUTION_METHOD);
    fclose(fileID);

else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end
