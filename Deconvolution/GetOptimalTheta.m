function theta = GetOptimalTheta(arr_fx,vDeg_arr_fx)

f = [1 -1 0];

nPolys_arr_fx = size(arr_fx,1);

mat_Ai = cell(nPolys_arr_fx,1);
mat_Bi = cell(nPolys_arr_fx,1);

for i = 0:1:nPolys_arr_fx-1

    % Get degree of the ith polynomial f_{i}(x)
    m = vDeg_arr_fx(i+1);
    
    mat_Ai{i+1,1} = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)'];

    mat_Bi{i+1,1} = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)'];
    
    vLambda{i+1,1} = log10(abs(arr_fx{i+1}));
    
end

Part_A = cell2mat(mat_Ai);
Part_B = cell2mat(mat_Bi);

A = [Part_A ; Part_B]; 

b = [cell2mat(vLambda) ; -1.*cell2mat(vLambda)];


x = linprog(f,-A,-b);

try
    theta = 10^x(3);

catch
    fprintf('Failed to optimize\n')
    theta = 1;

end

end