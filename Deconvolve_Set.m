function [h] = Deconvolve_Set(arr_fx)
% Deconvolve the set of polynomials q{i}
% 

% Get the number of polynomials in q
[~,nPolys] = size(arr_fx);

global SETTINGS

switch SETTINGS.DECONVOLVE_METHOD
    case 'Separate'
        for i = 1:1:nPolys-1
            h{i} = Deconvolve(arr_fx{i},arr_fx{i+1}) ;
        end
        
    case 'Batch'
        h = Deconvolve_Batch(arr_fx);
        

    otherwise
        error('SETTINGS.DECONVOLVE_METHOD is either Separate or Batch')
end
end


