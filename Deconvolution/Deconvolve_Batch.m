function arr_hx = Deconvolve_Batch(arr_fx)
% DECONVOLVE_BATCH Given the set of polynomials f_{0},...,f_{1}. compute
% the series of deconvolutions h_{1} = f_{1}/f_{0} h_{2} = f_{2}/f_{1},...
% Perform the deconvolutions by using the structure
% diag [ C(f_{1}) C(f_{2}) ... ] [h1 h2 ...]^{T} = [f_{0} f_{1} ...]^{T}
%
% % Inputs
%
% arr_fx : Array of polynomials f_{i}(x)
%
% % Outputs.
%
% arr_hx : Array of polynomials h_{i}(x)


% Get the number of polynomials in the set arr_fx
nPolys_arr_fx = size(arr_fx,1);
nPolys_arr_hx = size(arr_fx,1) - 1;

% Get the degree m_{i} of each of the polynomials f_{i} and store in a
% vector.
vDeg_arr_fx = zeros(nPolys_arr_fx,1);
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% Get the degrees n{i} of polynomials h_{i} = f_{i}/f_{i+1}.
vDeg_arr_hx = zeros(nPolys_arr_fx-1,1);
for i = 1:1:nPolys_arr_hx
    vDeg_arr_hx(i) = vDeg_arr_fx(i)-vDeg_arr_fx(i+1);
end

% Obtain theta such that the ratio of max element to min element is
% minimised

%
% y - Preprocess
% n - Dont preprocess
global SETTINGS
SETTINGS.PREPROC_DECONVOLUTIONS;

if( SETTINGS.PREPROC_DECONVOLUTIONS)
    
    theta = GetOptimalTheta(arr_fx,vDeg_arr_fx);
else
    theta = 1;
end


% Initialise a cell-array for f(w)
arr_fw = cell(size(arr_fx,1),1);

% for each f_{i} get fw_{i}
for i = 1:1:size(arr_fx,1)
    arr_fw{i} = GetWithThetas(arr_fx{i},theta);
end

RHS_vec = real(BuildRHSF(arr_fw));
DCQ = BuildC(arr_fw);

% Solve h_{0} for initial values of h
hw_vec = SolveAx_b(DCQ,RHS_vec);
v_h = hw_vec;

% Split vec h in to an array of polynomials.
arr_hw = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_hx
    
    % Get degree of h{i}
    deg_hw = vDeg_arr_hx(i);
    
    % Get coefficients of h_{i} from the solution vector
    arr_hw{i} = hw_vec(1:deg_hw+1);
    
    % Remove the coefficients from the solution vector
    hw_vec(1:deg_hw+1) = [];
end

% Remove thetas
arr_hx = cell(nPolys_arr_hx,1);
for i = 1:1: nPolys_arr_fx-1;
    arr_hx{i} = GetWithoutThetas(arr_hw{i},theta);
end

end


function f = BuildRHSF(fw_array)
% Build the vector f such that it contains the elements of
% Rhs f = [f_{0},...,f_{n-1}]
%
%
% fw = array of vectors f_{0},...,f_{n}
%

% Initialise empty vector.
f = [];

% for each vector f f_{0},...,f_{n-1} in fw_array, add to right hand
% side vector

for i=1:1:length(fw_array)-1
    f = [f;fw_array{i}];
end

end


function C = BuildC(arr_fx)
% Build the matrix C for the series of polynomial deconvolution
% f1/f2 = h1, f2/f3 = h2, ... , where the solutions h_{i} are given in the
% vector x of the expression Cx = b.
% The matrix  is a strucutred matrix of coefficients of polynomials
% f_{2},...,f_{n} where n is the number of polynomials in the set f_{i}. C
% has a block diagonal structure, where each block T_{i} on the diagonal is a
% toeplitz matrix of coefficients of f_{i+1}
%
% % Inputs.
%
% arr_fx = array of polynomials f(x)
%
% % Outputs.
%
% C : Matrix DCQ containing coefficients of f_{i} as described above.

% Get the number of polynomials in the array f_{i}
nPolys_f = length(arr_fx);

% Initialise an array of cells to store the submatrices which will make up
% DCQ
T = cell(nPolys_f-1,1);

% For each polynomial f_{i}
for i = 1:1:nPolys_f-1
    
    % Get the polynomial f_{i} = set_f{i+1}
    fw = arr_fx{i+1};
    
    % Get the polynomial f_{i-1} = set_f{i}
    fw_prev = arr_fx{i};
    
    % Get degree of polynomial f_{i} = m_{i}
    deg_fw = GetDegree(fw);
    
    % Get the degree of polynomial f_{i-1}
    deg_fw_prev = GetDegree(fw_prev);
    
    % Get the degree of the polynomial h_{i}
    deg_hw = deg_fw_prev - deg_fw;
    
    % Build the Matrix T(f)
    T{i} = BuildT1(fw,deg_hw);
end

% Build the Coefficient Matrix C of all matrices c
C = blkdiag(T{1:length(T)});

end

function Y = BuildY(set_hw,m)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

for i = 1:1:length(set_hw)
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h(w)
    hw = set_hw{i};
    
    % Get degree of f_{i}
    deg_fw = m(i+1);
    
    y{i} = real(BuildY1(hw,deg_fw));
end

%Build the Coefficient Matrix C
num_Rows = 0;
for i = 1:length(m)-1
    num_Rows = num_Rows + 1 + (m(i));
end
cols = (m(1)+1);

xx = zeros(num_Rows,cols);
Y = blkdiag( y{1:length(y)});
Y = [xx Y];

Y = Y;


end


function Y1 = BuildY1(hx,m1)
% Construct a partition Y1 of the matrix Y.
Y1 = BuildT1(hx,m1);
end
