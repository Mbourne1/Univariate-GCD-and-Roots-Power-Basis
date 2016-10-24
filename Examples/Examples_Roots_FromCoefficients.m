
function fx = Examples_Roots_FromCoefficients(ex_num)

addpath('../Examples');

% Get the symbolic factors of f(x) and corresponding multiplicities
f_root_sym_mult_array = Univariate_Roots_Examples(ex_num);

% Get the coefficients of f(x) as a vector.
fx = GetCoefficientsFromSymbolicRoots(f_root_sym_mult_array);


end

