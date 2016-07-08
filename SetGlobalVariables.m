function [] = SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method, bool_alpha_theta, low_rank_approx_method)
global SETTINGS

SETTINGS.PROBLEM_TYPE = problem_type;
SETTINGS.EX_NUM = ex_num;
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;
SETTINGS.MEAN_METHOD = mean_method;
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;
SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;
SETTINGS.SEED = 1024;

SETTINGS.THRESHOLD = 0.5;
SETTINGS.THRESHOLD_RANK = 1e-5;

SETTINGS.PLOT_GRAPHS = 'n';

SETTINGS.MAX_ERROR_SNTLN = 1e-16;
SETTINGS.MAX_ITE_SNTLN = 100;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-10;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;

%
% Separate
% Batch
%
SETTINGS.DECONVOLUTION_METHOD = 'Separate';

%
% From ux
% From Deconvolutions
%
SETTINGS.ROOTS_UX = 'From Deconvolutions';
end




