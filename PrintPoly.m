function [] = PrintPoly(fx,f)

[r,~] = size(fx);
m = r-1;

str = sprintf('%s(x) = ',f);


for i = 0:1:m

    % if the first coefficient, then dont preceed with plus
    if i == 0
        str_sign = '';
        str_monomial = '';
    else 
        % prefix with plus or minus depending on sign of coefficient
        if fx(i+1) >= 0
            str_sign = '+';
        else
            str_sign = '';
        end
        str_monomial = sprintf('x^{%i}',i);
    end
    
    
    
    str_temp = sprintf(' %2.4f ',fx(i+1));
    
    str = strcat(str,str_sign,str_temp,str_monomial);
    
    
end
   
fprintf(str)
fprintf('\n')
end