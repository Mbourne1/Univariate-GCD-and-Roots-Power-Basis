function [arr_hx] = Deconvolve_Set(arr_fx,DECONVOLUTION_METHOD)
% Deconvolve_Set
% Deconvolve the set of polynomials f_{i}(x), where the polynomails
% f_{i}(x) are outputs of a sequence of GCD computations in the Tobey and
% Horowitz algorithm.
%
% Inputs.
%
%
% arr_fx : Array of polynomials f_{i}(x)
%
% DECONVOLUTION_METHOD : (String)
%
% Outputs
%
%
% arr_hx : Array of polynomials h_{i}(x)

% Get the number of polynomials in array of f_{i}
nPolys_fx = size(arr_fx,1);

switch DECONVOLUTION_METHOD
    
    case 'Separate'
        
        arr_hx = cell(1,nPolys_fx-1);
        
        for i = 1:1:nPolys_fx-1
            arr_hx{i} = Deconvolve(arr_fx{i},arr_fx{i+1}) ;
        end
        
    case 'Batch'
        
        arr_hx = Deconvolve_Batch(arr_fx);
        
    case 'Batch With STLN'
        
        error([mfilename ' : Not yet completed Batch with STLN'])
        arr_hx = Deconvolve_Batch(arr_fx);
        
    case 'Batch Constrained'
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys_fx,1);
        for i = 1:1:nPolys_fx
            vDegt_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDeg_hx = diff(vDegt_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDeg_wx = diff([vDeg_hx; 0]);
        
        
        vMult = find(vDeg_wx~=0);
        
        arr_hx = Deconvolve_Batch_Constrained(arr_fx,vMult);
        
    case 'Batch Constrained With STLN'
        
        
        error([mfilename ' : ' 'Not Yet completed Batch Constrained with STLN']);
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys_fx,1);
        for i = 1:1:nPolys_fx
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


