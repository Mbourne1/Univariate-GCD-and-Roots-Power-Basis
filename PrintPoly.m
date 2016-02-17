function [] = PrintPoly(fx,f)

[r,~] = size(fx);
m = r-1;

str = sprintf('%s(x) = ',f);


for i = 0:1:m

    
    str_temp = sprintf(' + %2.4f x^{%i}  ',fx(i+1),i);
    
    str = strcat(str,str_temp);
    
    
end
   
fprintf(str)
fprintf('\n')
end