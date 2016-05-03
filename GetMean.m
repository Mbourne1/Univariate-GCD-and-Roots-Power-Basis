function lambda = GetMean(fx,n_k)
% Compute the mean of the non-zero entries of C_{n-k}(f).


global MEAN_METHOD
switch MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        % Get geometric mean
        lambda  = geomean(abs(fx(fx~=0)));       
    otherwise
        error('err')
        
end

end