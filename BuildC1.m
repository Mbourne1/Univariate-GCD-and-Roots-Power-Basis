function [C_f] = BuildC1(fx,n,k)
% Given the polynomial f(x) build the partition of the Sylvester matrix.
%
% Inputs
%
% fx :
%
% n  :
%
% k  :
%

%%

% Get degree of polynomial f
[r,~] = size(fx);
m = r -1;

C_f = zeros(m+n-k+1,n-k+1);

% for each column
for i = 0:1:n-k
   C_f(:,i+1) = [...
       zeros(i,1);
       fx;
       zeros((m+n-k+1)-(m+1)-i,1);...
       ];
end

end