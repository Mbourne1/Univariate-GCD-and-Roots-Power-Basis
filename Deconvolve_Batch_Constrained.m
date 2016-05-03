function [h] = Deconvolve_Batch_Constrained(arr_fx,vMult)

% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f(x). vMult = [m_{1}, m_{2} ,..., m_{n}]

% The division f_{0}/f_{1},...,f_{m_{1}}-1/f_{m_{1}} all have the same solution

% Get the set of polynomials f_{0} to f_{m_{1}-1}
for i = 1:1:length(vMult)
    if (i == 1)
        m1 = 1; 
    else 
        m1 = vMult(i-1);
    end
    
    m2 = vMult(i);
    subset_f{i} = arr_fx(1,m1:m2)
end


for i = 1:1:length(vMult)
    
    if i > 1
        old_mult = vMult(i-1)
    else
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    Cf{i} = [];
    
    fprintf('all polys from f(%i) to f(%i)',old_mult+1,new_mult+1)
    
    
    for j = old_mult+1+1:1:new_mult+1
        
        % Get the degree of the previous f_{i}(x)
        fx_prev = arr_fx{j-1};
        deg_fx_prev = GetDegree(fx_prev);
        
        % Get the degree of the polynomial f_{i}(x)
        fx = arr_fx{j}
        deg_fx = GetDegree(fx);
        
        
        Tf{j} = BuildT1(fx, deg_fx_prev - deg_fx );
        Cf{i} = [Cf{i} ; Tf{j}];
    end
    
    display(Cf)
    % Stack all T(f_{j})
    
    
end

LHS_Matrix = blkdiag(Cf{:});

RHS_vec = [];
% % Build the RHS vector
for i = 1:1:length(arr_fx)-1
    RHS_vec = [RHS_vec ; arr_fx{i}];
end


x = SolveAx_b(LHS_Matrix,RHS_vec);

x_temp = x;
for i = 1:1:length(vMult)
    mult = vMult(i)
    deg = GetDegree(arr_fx{mult}) - GetDegree(arr_fx{mult+1})
    h{i} = x_temp(1:deg+1)
    x_temp(1:deg+1) = [];
end


% Get Residual
residual = RHS_vec - (LHS_Matrix*x)
display(norm(residual))


end