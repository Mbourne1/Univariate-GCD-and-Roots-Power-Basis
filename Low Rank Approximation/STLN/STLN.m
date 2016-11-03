function [fx_lr,gx_lr] = STLN(fx,gx,k)
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
% fx_lr : Coefficients of f(x) after refinement f(x) + \delta
%
% gx_lr : Coefficients of g(x) after refinement g(x) + \delta

global SETTINGS

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get the derivative of f(x)
n = GetDegree(gx);

% Initialise the vector of perturbations zf(x)
zf = zeros(m+1,1);

% Initialise the vector of perturbations zg(x)
zg = zeros(n+1,1);

z = [zf ; zg];

% Build the k'th subresultant
T1_fx = BuildT1(fx,n-k);
T1_gx = BuildT1(gx,m-k);
Sk_fg = [T1_fx T1_gx];

% Build the matrix E_{t}(z)
T1_zf = BuildT1(zf,n-k);
T2_zg = BuildT1(zg,m-k);
Sk_zfzg = [T1_zf T2_zg];

% Get the index of the optimal colummn for removal
[~,idx_col] = GetMinDistance(Sk_fg);

% Get A_{t} the LHS matrix, equivalent to S_{t} with the optimal column
% removed.
At_fg = Sk_fg;
At_fg(:,idx_col) = [];

% Get c_{t} the removed column of S_{t} to form A_{t}.
ck = Sk_fg(:,idx_col);

% Get E_{t}, the matrix of structured perturbations corresponding to A_{t}.
Ak_zfzg = Sk_zfzg;
Ak_zfzg(:,idx_col) = [];

% Get h_{t}, the vector of structured perturbations corresponding to c_{t}
hk = Sk_zfzg(:,idx_col);

% Build Pt
Pk = BuildP(idx_col,m,n,k);

% Get the solution vector x of A_{t}x = c_{t}.
x_ls = SolveAx_b(At_fg,ck);

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ck + hk) - At_fg*x_ls;

% Build the matrix Y_{t}
vec_x = [x_ls(1:idx_col-1) ; 0 ; x_ls(idx_col:end)];

Yk = BuildY(vec_x,m,n,k);

% Build the matrix C for LSE Problem
H_z = Yk - Pk;
H_x = At_fg + Ak_zfzg;

C = [H_z H_x];

% Build the matrix E for LSE Problem
E = eye(2*m+2*n-2*k+3);


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


% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)./ norm(ck);



while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITE_SNTLN
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE(E,f,C,res_vec);
        
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    delta_zk        = y_lse(1:m+n+2,1);
    
    % Update z and x
    z = z + delta_zk;
    
    % Split z into z_f and z_g
    zf = z(1:m+1);
    zg = z(m+2:end);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    Sk_zfzg = BuildT(zf,zg,k);
    
    % Get the matrix E_{t} with optimal column removed
    Ak_zfzg = Sk_zfzg;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    hk = Sk_zfzg(:,idx_col);
    
   
    delta_xk        = y_lse((m+n+3):(2*m+2*n-2*k+3),1);
    x_ls = x_ls + delta_xk;
    vec_x = [x_ls(1:idx_col-1) ; 0 ; x_ls(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yk = BuildY(vec_x,m,n,k);
    
    % Get the residual vector
    res_vec = (ck+hk) - ((At_fg+Ak_zfzg)*x_ls);
    
    % Update the matrix C
    H_z = Yk - Pk;
    H_x = At_fg + Ak_zfzg;
    C = [H_z H_x];
    
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ck + hk) ;
    
end


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
    fx_lr = fx + zf;
    gx_lr = gx + zg;
else
    % Do nothing
    fx_lr = fx;
    gx_lr = gx;
end

fprintf([mfilename ' : ' sprintf('Required number of iterations : %i \n',ite)])



end


