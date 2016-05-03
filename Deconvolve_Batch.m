function h = Deconvolve_Batch(set_f)
% Given the set of polynomials f_{0},...,f_{1}. compute the series of
% deconvolutions h_{1} = f_{1}/f_{0} h_{2} = f_{2}/f_{1},...
% Perform the deconvolutions by producing the structure 
% diag [ C(f_{1}) C(f_{2}) ... ] [h1 h2 ...]^{T} = [f_{0} f_{1} ...]

nPolys_f = length(set_f);

% Get the degree of each of the polynomials
m = zeros(nPolys_f,1);
for i = 1:1:nPolys_f
    m(i) = GetDegree(set_f{i});
end


% Get values of n{i}: The degree of polynomials h_{i} for i = 1....d
n = zeros(1,nPolys_f-1);
for i = 1:1:nPolys_f-1
    n(i) = m(i)-m(i+1);
end

% for each polynomial f_{i}

T1 = cell(nPolys_f-1,1);

for i = 1:1:nPolys_f-1
    
    % Get the polynomial f_{i} = set_f{i+1} 
    fw = set_f{i+1};
    
    % Get the polynomial f_{i-1} = set_f{i}
    fw_prev = set_f{i};
    
    % Get degree of polynomial f_{i} = m_{i}
    deg_fw = GetDegree(fw);
    
    % Get the degree of polynomial f_{i-1}
    deg_fw_prev = GetDegree(fw_prev);
    
    % Get the degree of the polynomial h_{i}
    deg_hw = deg_fw_prev - deg_fw;
    
    % Build the Matrix T(f)
    T1{i} = BuildT1(fw,deg_hw);
end

T = blkdiag(T1{1:length(T1)});

RHS_vec = [];
% Build the vector [f_{0},...,f_{}-1]
for i = 1:1:nPolys_f - 1;
    % Get the polynomial f
    fw = set_f{i};
    
    RHS_vec = [RHS_vec; fw];
end

vec_h = SolveAx_b(T,RHS_vec);

% Split vec h in to an array of polynomials.
for i = 1:1:nPolys_f-1
    
    % Get degree of h{i}
    deg_hw = n(i);
    
    % Get coefficients of h_{i} from the solution vector
    h{i} = vec_h(1:deg_hw+1);
    
    % Remove the coefficients from the solution vector
    vec_h(1:deg_hw+1) = [];
end

end

