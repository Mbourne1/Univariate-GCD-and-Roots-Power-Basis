function St = BuildSubresultant(fw,gw,t,alpha)
% Build the Sylvester Subresultant matrix S_{t}(f,g)

% Get degree of polynomial f(w)
[rows_f, ~] = size(fw);
m = rows_f - 1 ;

% Get degree of polynomial g(w)
[rows_g, ~] = size(gw);
n = rows_g - 1 ;

C1 = BuildC1(fw,n-t);
C2 = BuildC1(gw,m-t);


St = [C1 alpha.*C2];

end