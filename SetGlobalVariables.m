function [] = SetGlobalVariables(mean_method, bool_alpha_theta, low_rank_approx_method)

global MEAN_METHOD
MEAN_METHOD = mean_method;

global BOOL_ALPHA_THETA
BOOL_ALPHA_THETA = bool_alpha_theta;

global LOW_RANK_APPROXIMATION_METHOD
LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

% Initialise the Seed for random noise generation
global SEED
SEED = 1024;

% Set the Threshold value for the GetGCDDegree function
global THRESHOLD
THRESHOLD = 3;

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

global BOOL_NOISE
BOOL_NOISE = 'y';

global MAX_ERROR_SNTLN
global MAX_ITE_SNTLN

MAX_ERROR_SNTLN = 1e-14;
MAX_ITE_SNTLN = 100;

global DECONVOLVE_METHOD 
DECONVOLVE_METHOD = 'Batch';


end