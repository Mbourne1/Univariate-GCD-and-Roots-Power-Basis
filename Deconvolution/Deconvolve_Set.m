function [arr_hx] = Deconvolve_Set(arr_fx,DECONVOLUTION_METHOD)
% Deconvolve the set of polynomials f_{i}(x)


% Get the number of polynomials in q
[~,nPolys] = size(arr_fx);

switch DECONVOLUTION_METHOD
    case 'Separate'
        
        arr_hx = cell(1,nPolys-1);
        
        for i = 1:1:nPolys-1
            arr_hx{i} = Deconvolve(arr_fx{i},arr_fx{i+1}) ;
        end
        
    case 'Batch'
        arr_hx = Deconvolve_Batch(arr_fx);
        
    case 'Batch With STLN'
        arr_hx = Deconvolve_Batch(arr_fx);
        
    case 'Batch Constrained'
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys,1);
        for i = 1:1:nPolys
            vDegt_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDeg_hx = diff(vDegt_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDeg_wx = diff([vDeg_hx; 0]);
        
        
        vMult = find(vDeg_wx~=0);
        
        arr_hx = Deconvolve_Batch_Constrained(arr_fx,vMult);
        
    case 'Batch Constrained With STLN'
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys,1);
        for i = 1:1:nPolys
            vDegt_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDeg_hx = diff(vDegt_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDeg_wx = diff([vDeg_hx; 0]);
        
        
        vMult = find(vDeg_wx~=0);
        
        arr_hx = Deconvolve_Batch_Constrained(arr_fx,vMult);
        
    otherwise
        error('SETTINGS.DECONVOLVE_METHOD is either Separate or Batch')
end
end


