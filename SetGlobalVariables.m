function [] = SetGlobalVariables(bool_preproc,low_rank_approx_method)

global BOOL_PREPROC
BOOL_PREPROC = bool_preproc;

global LOW_RANK_APPROXIMATION_METHOD
LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

% Initialise the Seed for random noise generation
global SEED
SEED = 1024;

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

global BOOL_NOISE
BOOL_NOISE = 'y';

global MAX_ERROR_SNTLN
global MAX_ITE_SNTLN

MAX_ERROR_SNTLN = 1e-14;
MAX_ITE_SNTLN = 100;

end