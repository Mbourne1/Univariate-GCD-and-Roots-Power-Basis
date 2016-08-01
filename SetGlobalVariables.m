function [] = SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method, bool_alpha_theta, low_rank_approx_method)
global SETTINGS

% % Problem Type
% Roots
% GCD
%
SETTINGS.PROBLEM_TYPE = problem_type;

% Example Number
SETTINGS.EX_NUM = ex_num;

% Minimum/Maximum noise level
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

% Method of computing mean
SETTINGS.MEAN_METHOD = mean_method;

% Include/Exclude preprocessing 
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;

% Low rank approximation method
%
% STLN
% SNTLN
SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

% Seed for noise generation
SETTINGS.SEED = 1024;

% Threshold defined for use in computing the degree of the GCD
SETTINGS.THRESHOLD = 0.5;
SETTINGS.THRESHOLD_RANK = 1e-5;

% Include/Exclude plotting of graphs
SETTINGS.PLOT_GRAPHS = 'y';

% Settings for SNTLN/STLN
SETTINGS.MAX_ERROR_SNTLN = 1e-15;
SETTINGS.MAX_ITE_SNTLN = 100;

% Settings for deconvolution
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-15;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;

% (Roots Only)
% Make use of precalculated limits on the upper and lower bounds of the
% degree t of the GCD(f,f').
SETTINGS.BOOL_LIMITS = 'y';

% Deconvolution Method
% Method for performing deconvolutions.
%
% Batch
% Separate
% Batch With STLN
% Batch Constrained
% Batch Constrained With STLN
%
SETTINGS.DECONVOLUTION_METHOD_FX_HX = 'Separate';

SETTINGS.DECONVOLUTION_METHOD_HX_WX = 'Separate';

%
% From ux
% From Deconvolutions
%
SETTINGS.ROOTS_UX = 'From Deconvolutions';
end




