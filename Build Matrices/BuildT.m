function Sk = BuildT(fx, gx, k)
% Build the Sylvester Subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% k : Index of Sylvester Subresultant matrix to be constructed.


% Get degree of polynomial f(w)
m = GetDegree(fx);

% Get degree of polynomial g(w)
n = GetDegree(gx);

% Build the matrix C(f)
C1 = BuildT1(fx,n-k);

% Buiild the matrix C(g)
C2 = BuildT1(gx,m-k);


Sk = [C1 C2];

end