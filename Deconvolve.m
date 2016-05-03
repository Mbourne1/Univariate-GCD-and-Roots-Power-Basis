function [h] = Deconvolve(arr_fx)
% Deconvolve the set of polynomials q{i}
% 

% Get the number of polynomials in q
[~,nPolys] = size(arr_fx);

global DECONVOLVE_METHOD

switch DECONVOLVE_METHOD
    case 'Single'
        for i = 1:1:nPolys-1
            h{i} = Deconvolve_Separate(arr_fx{i},arr_fx{i+1}) ;
        end
        
    case 'Batch'
        h = Deconvolve_Batch(arr_fx);
        

    otherwise
        error('err')
end
end

function [h] = Deconvolve_Separate(f,g)

% Get degree of polynomial f(x).
m = GetDegree(f);

% Get degree of polynomial g(x).
n = GetDegree(g);

% Build the matrix C(g
C_g = BuildT1(g,m-n);

% Solve C(g)*h = f
h = SolveAx_b(C_g,f);

end
