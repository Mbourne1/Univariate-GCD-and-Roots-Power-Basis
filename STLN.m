function [fx,gx] = STLN(fx,gx,t)
% Perform Structured Total Least Norm to obtain a low rank approximation 
% of the t-th Sylvester matrix. Note this is a linear problem, any
% alpha and theta values are already included in f(x) and g(x).
% 
%
% % Inputs.
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
%
% t : Degree of GCD(f(x),g(x))
%
% % Outputs.
%
% fx : Coefficients of f(x) after refinement f(x) + \delta
%
% gx : coefficients of g(x) after refinement g(x) + \delta


% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get the derivative of f(x)
n = GetDegree(gx);

% Initialise the vector of perturbations zf(x)
zf = zeros(m+1,1);

% Initialise the vector of perturbations zg(x)
zg = zeros(n+1,1);

z = [zf ; zg];

% Build the t'th subresultant
C1 = BuildT1(fx,n-t);
C2 = BuildT1(gx,m-t);
St = [C1 C2];

% Build the matrix E_{t}(z)
B1 = BuildT1(zf,n-t);
B2 = BuildT1(zg,m-t);
Bt = [B1 B2];

% Get the index of the optimal colummn for removal
[~,colIndex] = GetMinDistance(St);

% Get A_{t} the LHS matrix, equivalent to S_{t} with the optimal column
% removed.
At = St;
At(:,colIndex) = [];

% Get c_{t} the removed column of S_{t} to form A_{t}.
ct = St(:,colIndex);

% Get E_{t}, the matrix of strucured perturbations corresponding to A_{t}.
Et = Bt;
Et(:,colIndex) = [];

% Get h_{t}, the vector of strucutred perturbations corresponding to c_{t}
ht = Bt(:,colIndex);

% Build Pt
Pt = BuildPt(colIndex,m,n,t);

% Get the solution vector x of A_{t}x = c_{t}.
x_ls = SolveAx_b(At,ct);

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
g = (ct + ht) - At*x_ls;

% Build the matrix Y_{t}
x = [x_ls(1:colIndex-1) ; 0 ; x_ls(colIndex:end)];

Yt = BuildYt(x,m,n,t);

% Build the matrix C for LSE Problem
H_z = Yt - Pt;
H_x = At+Et;

C = [H_z H_x];

% Build the matrix E for LSE Problem
E = eye(2*m+2*n-2*t+3);
%E = blkdiag(eye(m+n+2),zeros(m+n-2*t+1,m+n-2*t+1))
% Build the matrix D which accounts for repetitions of z_{i} in B_{k}
%E = blkdiag(eye(n-t+1),eye(m-t+1));


% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    z;...
    x_ls;
    ];

% Set yy to be the vector which stores all cummulative perturbations.
yy = start_point;

% Set the initial value of vector p to be zero
f = -(yy-start_point);
%f = start_point

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(g)./ norm(ct);
global SETTINGS
while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITE_SNTLN

    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE(E,f,C,g);
    
    %y_lse = lsqlin(E,f,C,g);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    delta_zk        = y_lse(1:m+n+2,1);
    delta_xk        = y_lse((m+n+3):(2*m+2*n-2*t+3),1);
    
    % Update z and x
    z = z + delta_zk;
    x_ls = x_ls + delta_xk;
    
    % Split z into z_f and z_g
    zf = z(1:m+1); 
    zg = z(m+2:end);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    Bt = BuildT(zf,zg,t);
    
    % Get the matrix E_{t} with optimal column removed 
    Et = Bt;
    Et(:,colIndex) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    ht = Bt(:,colIndex);
    
    % Get the updated vector x
    x = ...
        [
        x_ls(1:colIndex-1); 
        0 ;
        x_ls(colIndex:end);
        ];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yt = BuildYt(x,m,n,t);
    
    % Get the residual vector
    g = (ct+ht) - ((At+Et)*x_ls);
    
    % Update the matrix C
    H_z = Yt - Pt;
    H_x = At + Et;
    C = [H_z H_x];
    
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(g)./norm(ct+ht) ;
    
end

global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s - Residuals',mfilename);
        figure('name',figure_name)
        hold on
        plot(log10(condition),'-s');
        hold off
        
    case 'n'
end

% If the final condition is less than the original, output the new values,
% otherwise output the old values for f(x) and g(x).
if (condition(ite) < condition(1))
    fx = fx + zf;
    gx = gx + zg;
else
  % Do nothing
end

fprintf([mfilename ' : ' sprintf('Required number of iterations : %i \n',ite)])

end


function Pt = BuildPt(idx_Col,m,n,t)
% Build the matrix P_{t}, where h_{t} = P_{t}z

if idx_Col <= n-t+1
    % First Partition
    i = idx_Col;
    Pt = ...
        [
            zeros(i-1,m+1)      zeros(i-1,n+1);
            eye(m+1,m+1)        zeros(m+1,n+1);
            zeros(n-t-i+1,m+1)  zeros(n-t-i+1,n+1);
        ];
else
    % Second Partition
    i = idx_Col - (n-t+1);
    Pt = ...
        [
            zeros(i-1,m+1)      zeros(i-1,n+1);
            zeros(n+1,m+1)      eye(n+1,n+1) 
            zeros(m-t-i+1,m+1)  zeros(m-t-i+1,n+1);
        ];
    
end

end

function Yt = BuildYt(x,m,n,t)
% Build the matrix Y_{t}, where Y_{t}z = E_{t}x

xa = x(1:n-t+1);
xb = x(n-t+2:end);

Y1 = BuildT1(xa,m);
Y2 = BuildT1(xb,n);

Yt = [Y1 Y2];


end