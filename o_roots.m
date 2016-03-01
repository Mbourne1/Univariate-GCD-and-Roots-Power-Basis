function [] = o_roots(ex_num,el,bool_sntln,seed)
% Calculate the roots of a polynomial f(x) by a series of methods.
%
%   Inputs.
%
%
%   ex_num :
%
%   el :
%
%   bool_sntln :

%%

% Global Variables
global BOOL_PREPROC
BOOL_PREPROC = 'n';

global BOOL_SNTLN
BOOL_SNTLN = bool_sntln;

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

global SEED
SEED = seed;

global MAX_ERROR_SNTLN
global MAX_ITE_SNTLN

MAX_ERROR_SNTLN = 1e-12;
MAX_ITE_SNTLN = 100;
%%
% Get polynomial f(x)

EXAMPLE_TYPE = 'FromRoots'

switch EXAMPLE_TYPE 
    case 'FromRoots'
        [       fx] = Examples_Roots(ex_num);
    case 'FromCoefficients'
        switch ex_num
            case '1'
                fx = [1; 1; -2; -2; 1; 1];
            case '2'
                fx = [-1; -1; 1; 1];
            case '3'
                fx = [-5; 1; -5; 1]
        end
    otherwise 
        error('Example polynomial is either from roots or from coefficients')
end


% Add Noise to f(x)
fx = Noise(fx,el);

%%
% Get roots by my method
o_roots_mymethod(fx);


% Get roots by matlab method
roots_calc = [roots(flipud(fx)) ones(length(roots(flipud(fx))),1)];
% Print the calculated roots and the corresponding multiplicities.
PrintRoots(roots_calc,'MATLAB METHOD');

% Get roots by zheng method
o_roots_multroot(flipud(fx));

end

