function arr_hx = Deconvolve_Batch_With_STLN(set_f)
% Given the set of polynomials f_{0},...,f_{1}. compute the series of
% deconvolutions h_{1} = f_{1}/f_{0} h_{2} = f_{2}/f_{1},...
% Perform the deconvolutions by producing the structure
% diag [ C(f_{1}) C(f_{2}) ... ] [h1 h2 ...]^{T} = [f_{0} f_{1} ...]

% Global Variables.
global SETTINGS


% Get the number of polynomials in the set set_f
nPolys_fx = length(set_f);

% let d be the number of deconvolutions = number of polys in set_f - 1;
d = nPolys_fx - 1;

% Get the degree m_{i} of each of the polynomials f_{i} and store in a
% vector.
m = zeros(nPolys_fx,1);
for i = 1:1:nPolys_fx
    m(i) = GetDegree(set_f{i});
end

% Get the degrees n{i} of polynomials h_{i} = f_{i}/f_{i+1}.
n = zeros(1,nPolys_fx-1);
for i = 1:1:nPolys_fx-1
    n(i) = m(i)-m(i+1);
end

% Define M to be the total number of all coefficients of the first d polynomials
% f_{0}...f_{d-1},
M = sum(m+1) - (m(end)+1);

% Define M1 to be the total number of all coefficients of polynomials
% f_{0},...,f_{d}
nCoefficients_fx = sum(m+1);

% Define N to be the number of coefficients of all h_{i}
nCoefficients_hx = sum(n+1);

% Obtain theta such that the ratio of max element to min element is
% minimised
%theta = getOptimalTheta(set_f,m);
theta = 1;

% Initialise a cell-array for f(w)
fw = cell(1,nPolys_fx);

% for each f_{i} get fw_{i}
for i = 1 : 1 : nPolys_fx 
    fw{i} = GetWithThetas(set_f{i},theta);
end

RHS_vec = real(BuildRHSF(fw));

Cf = BuildC(fw);

% Solve h_{0} for initial values of h
v_hx = SolveAx_b(Cf,RHS_vec);

% Split vec h in to an array of polynomials.
for i = 1:1:nPolys_fx-1
    
    % Get degree of h{i}
    deg_hw = n(i);
    
    % Get coefficients of h_{i} from the solution vector
    arr_hx{i} = v_hx(1:deg_hw+1);
    
    % Remove the coefficients from the solution vector
    v_hx(1:deg_hw+1) = [];
end

% Let z be  vectors of perturbations to polynomials fi such that
% z = [z{0} z{1} z{2} z{3} ... z{d}]
arr_z = cell(1,nPolys_fx);

for i = 1 : 1 : nPolys_fx
    arr_z{i} = zeros(1,m(i)+1);
end

% Build vector z, consisting of all vectors z_{i}
z_o = [arr_z{:}]';

% Build the Matrix P
P = [eye(M) zeros(M,nCoefficients_fx-M)];

% Build Matrix Y, where E(z)h = Y(h)z
Y = BuildY(arr_hx,m);

% Set the iteration counter.
ite = 1;


F = eye(nCoefficients_hx + nCoefficients_fx);

G = [Cf Y-P];

% Compute the first residual
res_vec = (RHS_vec + (P*z_o) - (Cf*v_hx));

zx_ite = z_o;

% Update Matrix Pz
Pz = P*zx_ite;

% Get the initial residual
condition(ite) = norm(res_vec)./norm(RHS_vec + Pz);

v_hx_ite = v_hx;

start_point = ...
    [
        v_hx;
        z_o;
    ];

s = yy - start_point;

yy = start_point;

% Perform iteration to obtain perturbations

while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(F,s,G,res_vec);
    
    yy = yy + y;
    
    % Output y gives delta h and delta z
    delta_h = y(1:nCoefficients_hx);
    delta_z = y(nCoefficients_hx+1:end);
    
    % Add structured perturbations to vector hx.
    v_hx_ite = v_hx_ite + delta_h;
    
    % Add structured perturbations to vector z.
    zx_ite = zx_ite + delta_z;
    
    % Seperate delta_z into its component vectors delta_z0 delta_z1,...,
    % delta_zd
    zz = zx_ite;
    
    arr_zx_ite = cell(1,nPolys_fx);
    for i = 1:1:nPolys_fx
        arr_zx_ite{i} = zz(1:m(i)+1);
        zz(1:m(i)+1) = [];
    end
    
    % Renew Matrix Pz
    Pz = P*zx_ite;
    
    % Increment s in LSE Problem
    s = -(yy - start_point);
    
    % Copy vector hx_ite 
    hx_temp = v_hx_ite;
    
    % Move individual vectors hi into variable size array, emptying
    % hh
    for i = 1:1:length(n)
       
        % Get set of coefficients of polynomial h_{i}(x)
        arr_hx{i} = hx_temp(1:n(i)+1);
        
        % Remove coeffficients from hx_temp
        hx_temp(1:n(i)+1) = [];
    end
    
    % Build iterative Y, where Y
    Y = BuildY(arr_hx,m);
    
    
    % Add the structured perturbations to improved fx array.
    arr_new_fx = cell(1,nPolys_fx);
    for i = 1:length(fw)
        arr_new_fx{i} = fw{i} + arr_zx_ite{i};
    end
    
    % Build the matrix CE = C(f) + E(z)
    CE = BuildC(arr_new_fx) ;
    
    % Build G
    G = [CE (Y-P)];
    
    % Calculate residual and increment t in LSE Problem
    res_vec = ((RHS_vec+Pz) - (CE*v_hx_ite));
    
    
    % Increment iteration number
    ite = ite + 1;
    
    condition(ite) = norm(res_vec)./norm(RHS_vec+Pz);
    
    
end



% Print outputs to command line
fprintf([mfilename ' : ' 'Performed Deconvolutions\n'])
fprintf([mfilename ' : ' sprintf('Iterations required for Batch Deconvolution %i\n', ite)])

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Condition',mfilename);
        figure('name',figure_name)
        hold on
        plot(log10(condition),'-s')
        hold off
    case 'n'
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

function Y = BuildY(arr_hx,m)
% Build the coefficient matrix Y. This is the change of variable such
% that
% E(z)*h = Y(h)*f

for i = 1:1:length(arr_hx)
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h_{i}(x)
    hx = arr_hx{i};
    
    % Get degree of f_{i}(x)
    deg_fx = m(i+1);
    
    Ch{i} = (BuildY1(hx,deg_fx));
end

% Build the Coefficient Matrix C

% Must include a diagonal section of zeros, since C*f =
num_Rows = 0;
for i = 1:length(m)-1
    num_Rows = num_Rows + 1 + (m(i));
end
cols = (m(1)+1);

xx = zeros(num_Rows,cols);
Y = blkdiag( Ch{1:length(Ch)});
Y = [xx Y];

Y = Y;


end


function Y1 = BuildY1(hx,m1)
% Construct a partition Y1 of the matrix Y.
Y1 = BuildT1(hx,m1);
end
