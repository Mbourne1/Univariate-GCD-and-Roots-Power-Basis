function [] = PrintCoefficientsBivariate(fxy,f)
% Given the polynomial f(x,y) print out the polynomial.

[r,c] = size(fxy);
m1 = r - 1 ;
m2 = c - 1 ;

str = sprintf('%s(x) = ',f)

% for each row
for i = 0:1:m1
    % for each column
    for j = 0:1:m2
        % 
        temp_str = sprintf('+ %2.4f x^{%i} y^{%i}',fxy(i+1,j+1),i,j);
        str = strcat(str,temp_str);
    end
end

fprintf(str)
fprintf('\n')
