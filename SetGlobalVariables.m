function [] = SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method, bool_alpha_theta, low_rank_approx_method)
global SETTINGS

% % Problem Type
% Roots
% GCD
%
SETTINGS.PROBLEM_TYPE = problem_type;

% Example Number
SETTINGS.EX_NUM = ex_num;

% Seed for noise generation
SETTINGS.SEED = 1024;

% Minimum/Maximum noise level
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

% Include/Exclude plotting of graphs
SETTINGS.PLOT_GRAPHS = 'y';


%--------------------------------------------------------------------------
%
%           SETTINGS : PREPROCESSING
%
%--------------------------------------------------------------------------



% Method of computing mean
SETTINGS.MEAN_METHOD = mean_method;

% Include/Exclude preprocessing 
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;



%--------------------------------------------------------------------------
%
%       SETTINGS : DEGREE COMPUTATION
%
%
%
%--------------------------------------------------------------------------



% Threshold defined for use in computing the degree of the GCD
SETTINGS.THRESHOLD = 3.1;
SETTINGS.THRESHOLD_RANK = 1e-5;

% (Roots Only)
% Make use of precalculated limits on the upper and lower bounds of the
% degree t of the GCD(f,f').
SETTINGS.BOOL_LIMITS = 'y';


%--------------------------------------------------------------------------
% 
%       SETTINGS : Low rank approximation method
%
%
%--------------------------------------------------------------------------

%
% STLN
% SNTLN
SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

% Settings for SNTLN/STLN
SETTINGS.MAX_ERROR_SNTLN = 1e-12;
SETTINGS.MAX_ITE_SNTLN = 100;





%-------------------------------------------------------------------------
%
% 
%               DECONVOLUTION SETTINGS
%
%-------------------------------------------------------------------------

% Deconvolution Method
% Method for performing deconvolutions.
%
% Batch
% Separate
% Batch With STLN
% Batch Constrained
% Batch Constrained With STLN
%


% Settings for deconvolution
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-13;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;


% Set deconvolution method for computing set of polynomials h_{i}(x) from
% the set of polynomials f_{i}(x)

%
% 'Separate'
% 'Batch'
% 'Batch With STLN'
% 'Batch Constrained'
% 'Batch Constrained With STLN'
%

SETTINGS.DECONVOLUTION_METHOD_FX_HX = 'Batch Constrained With STLN';

% Set the deconvolution method for computing the set of polynomials
% w_{i}(x) from the set of polynomials h_{i}(x)
SETTINGS.DECONVOLUTION_METHOD_HX_WX = 'Separate';

SETTINGS.PREPROC_DECONVOLUTIONS = 'y';

%
% From ux
% From Deconvolutions
%
SETTINGS.ROOTS_HX = 'From ux';
end




