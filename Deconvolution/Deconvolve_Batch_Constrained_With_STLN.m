function [arr_hx] = Deconvolve_Batch_Constrained_With_STLN(arr_fx,vMult)
% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f_{0}(x). vMult = [m_{1}, m_{2} ,..., m_{n}]
% The division f_{0}/f_{1},...,f_{m_{1}-1} / f_{m_{1}} all have the same solution
%
%
% % Inputs.
%
% arr_fx : Array of polynomials f(x)
%
% vMult : Multiplicities of the factors of f(x) in ascending order.
%
% % Outputs.
%
% arr_hx : Array of polynomials h(x) where h_{i}(x) = f_{i-1}(x) / f_{i}(x)


% Global Variables.
global SETTINGS


% Get the number of polynomials in the set set_f
nPolys_arr_fx = length(arr_fx);

% let d be the number of deconvolutions = number of polys in set_f - 1;
d = nPolys_arr_fx - 1;

% Get the degree m_{i} of each of the polynomials f_{i} and store in a
% vector.
m = zeros(nPolys_arr_fx,1);
for i = 1:1:nPolys_arr_fx
    m(i) = GetDegree(arr_fx{i});
end

% Get the degrees n{i} of polynomials h_{i}(x) = f_{i}(x)/f_{i+1}(x).
n = m(1:end-1) - m(2:end);

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

% %
% %
% Build LHS Matrix
LHS_Matrix = BuildC(arr_fx,vMult);

% %
% %
% Build the RHS vector

% RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
% degree of any root of f_{0}(x).
RHS_vec = [];
for i = 1 : 1 : nPolys_arr_fx - 1
    RHS_vec = [RHS_vec ; arr_fx{i}];
end

% Get vector of coefficients of h_{i}(x) for all i.
v_px = SolveAx_b(LHS_Matrix,RHS_vec);

nCoefficients_px = length(v_px);

% Get unique multiplities from the multiplicity vector
unique_vMult = unique(vMult);

% Initialise a temporary vector
x_temp = v_px;

nEntries_arr_px = length(unique_vMult);


vDeg_px = zeros(nEntries_arr_px , 1);

for i = 1:1: nEntries_arr_px
    mult = unique_vMult(i);
    deg = GetDegree(arr_fx{mult}) - GetDegree(arr_fx{mult+1});
    vDeg_px(i) = deg;
    arr_px{i} = x_temp(1:deg+1);
    x_temp(1:deg+1) = [];
end


% %
% %
% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).

count = 1;
for i = 1:1:nEntries_arr_px
    
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hx{count} = arr_px{i};
        count = count + 1;
    end
    
end


% Build the array of polynomials z(x) which are the structured
% perturbations of the array of polynomials f(x).

nPolys_arr_fx = size(arr_fx,1);
arr_zx = cell(1,nPolys_arr_fx);

for i = 1 : 1 : nPolys_arr_fx
    arr_zx{i} = zeros(1,m(i) +1);
end

% Build the vector zx consisting of all vectors in arr_zx
v_zx = [arr_zx{:}]';

% Build the matrix P
P = [eye(M) zeros(M,nCoefficients_fx-M)];

% Compute the first residual
res_vec = RHS_vec + (P*v_zx) - (LHS_Matrix * v_px);

% Update Matrix P*z
Pz = P*v_zx;

% Build Matrix Y, where E(z)h = Y(h)z
Y = BuildY(arr_hx,m);

%--------------
for i = 1 : 1 : nPolys_arr_fx
    vec_fx = [RHS_vec; arr_fx{i}];
end
test1 = Y*vec_fx;
test2 = LHS_Matrix*v_px;
test1./test2;
%--------------


% Set the iteration number
ite = 1;

F = eye(nCoefficients_px + nCoefficients_fx);

G = [LHS_Matrix Y-P];

% Get the intial residual
condition(ite) = norm(res_vec)./norm(RHS_vec + Pz);

start_point = ...
    [
    v_px;
    v_zx;
    ];

yy = start_point;

s = -(yy - start_point);

while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS) && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    y = LSE(F,s,G,res_vec);
    
    yy = yy + y;
    
    delta_px = y(1:nCoefficients_px);
    delta_zx = y(nCoefficients_px+1:end);
    
    % Add structured perturbations to zx.
    v_zx = v_zx + delta_zx;
    
    %
    zz = v_zx;
    
    % Separate delta_z into its component polynomials zx{i} corresponding
    % to fx{i}.
    for i = 1:1:nPolys_arr_fx
        arr_zx{i} = zz(1:m(i)+1);
    end
    
    % Add structured perturbations to improved f
    for i = 1:1:length(arr_fx)
        arr_fx{i} = arr_fx{i} + arr_zx{i};
    end
    
    
    % RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
    % degree of any root of f_{0}(x).
    RHS_vec = [];
    for i = 1 : 1 : nPolys_arr_fx - 1
        RHS_vec = [RHS_vec ; arr_fx{i}];
    end

    % Get vector of coefficients of h_{i}(x) for all i.
    v_px = SolveAx_b(LHS_Matrix,RHS_vec);
    v_px = v_px + delta_px;
    
    % Build the matrix CE = C(f) + E(z)
    LHS_Matrix = BuildC(arr_fx,vMult);
    
    % Renew the matrix P*z
    Pz = P* v_zx;
    
    % Increment s in LSE Problem.
    s = -(yy-start_point);
    
    % Copy vector px_ite
    px_temp = v_px;
    
    % Move individual vectors p_{i}(x) into array of polynomials arr_px
    for i = 1:1:nEntries_arr_px
        
        % Get set of coefficients of polynomial h_{i}(x)
        deg = vDeg_px(i);
        arr_px{i} = px_temp(1:deg+1);
        
        % Remove coefficients from hx_temp
        px_temp(1:deg+1) = [];
    end
    
    % %
    % Copy px to obtain the set of polynomials hx
    count = 1;
    for i = 1:1:nEntries_arr_px
        
        if i == 1
            nReps = unique_vMult(i);
        else
            nReps = (unique_vMult(i) - unique_vMult(i-1));
        end
        
        for j = 1:1:nReps
            arr_hx{count} = arr_px{i};
            count = count + 1;
        end
        
    end
    
    % Build iterative Y
    Y = BuildY(arr_hx,m);
    
    % Build G
    G = [LHS_Matrix (Y-P)];
    
    % Calculate residual and increment t in LSE Problem
    res_vec = ((RHS_vec+Pz) - (LHS_Matrix*v_px));
    
    % Increment iteration number
    ite = ite + 1;
    
    condition(ite) = norm(res_vec)./norm(RHS_vec + Pz);
    
    
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
    
    % Get degree of f_{i}
    deg_fw = m(i+1);
    
    y{i} = real(BuildY1(hx,deg_fw));
end

% Build the Coefficient Matrix C

% Get number of columns for zero segment
nCols = (m(1)+1);


Y = blkdiag( y{1:length(y)});
nRows = size(Y,1);

Y = [zeros(nRows,nCols) Y];

Y = Y;


end

function Y1 = BuildY1(hx,m1)
% Construct a partition Y1 of the matrix Y.
Y1 = BuildT1(hx,m1);
end

function LHS_Matrix = BuildC(arr_fx,vMult)
% %
% %
% Build the LHS Coefficient matrix

% For each distinct hx, build the partition of the Matrix C(f_{m_{i}+1},...,f_{m_{i+1}})
% C(f_{1},...f_{m1}), C(f_{m1+1},...,f_{m2}),...
% C(f_{1},...,f_{m1}) = [T(f1) ; T(f2) ; ... T(fm1)]

% Get number of distinct polynomials in h_{i}
nDistinct_hx = length(vMult);

for i = 1:1:nDistinct_hx
    
    if i > 1
        old_mult = vMult(i-1);
    else % No previous multiplicity
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    Cf{i} = [];
    
    % for each polynomial f_{i} in the interval f_{m_{i-1}+1}...f_{m_{i}}
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get the degree of the previous polynomial f_{i-1}(x)
        fx_prev = arr_fx{j-1};
        deg_fx_prev = GetDegree(fx_prev);
        
        % Get the degree of the current polynomial f_{i}(x)
        fx = arr_fx{j};
        deg_fx = GetDegree(fx);
        
        % Get the degree of the polynomial h_{i} = f_{i-1}/f_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the Cauchy like matrix T_{m_{i} - m_{i-1}}(f_{i})
        Tf{j} = BuildT1(fx, deg_hx);
        
        % Stack beneath all other T_{f} which are multiplied by [_{i}(x)
        Cf{i} = [Cf{i} ; Tf{j}];
    end
    
    
end

LHS_Matrix = blkdiag(Cf{:});


end
