function St = BuildT(fx,gx,t)
% Build the Sylvester Subresultant matrix S_{t}(f,g)

% Get degree of polynomial f(w)
m = GetDegree(fx);

% Get degree of polynomial g(w)
n = GetDegree(gx);

% Build the matrix C(f)
C1 = BuildT1(fx,n-t);

% Buiild the matrix C(g)
C2 = BuildT1(gx,m-t);


St = [C1 C2];

end