

p = [1 5.503 9.765 7.647 2.762 0.37725];
q = [1 -2.993 -0.7745 2.0070 0.7605];

m = length(p) -1;
n = length(q) -1;

S = zeros(m+n,m+n);

for i = 0:1:n
   S(i+1,i+1:i+m+1) = p; 
end

for i = n+1:1:m+n+1
   S(i+1,i-n:i) = q; 
end

vSingularValues = svd(S)
