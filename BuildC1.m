function [C_f] = BuildC1(fx,n_k)
% Given the polynomial f(x) build the partition C1 of the Sylvester matrix.
%
% Inputs
%
% fx    :
%
% n_k   :
%

%%

% Get degree of polynomial f
[r,~] = size(fx);
m = r -1;

C_f = zeros(m+n_k+1,n_k+1);

% for each column
for i = 0:1:n_k
   C_f(i+1:i+m+1,i+1) = fx;
end

end