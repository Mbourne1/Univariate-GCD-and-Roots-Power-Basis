function [lambda,mu,alpha,theta] = Preproccess(fx,gx)
% Preprocess(fx,gx)
%
% Preprocess the input polynomials to obtain geometric means of the
% coefficients of f(x) and g(x), and optimal values of alpha and theta such
% that the ratio of entries in sylvester matrix S(f,g) is minimized.
%
% Inputs.
%
% fx : Vector of coefficients of polynomial f(x)
%
% gx : Vector of coefficients of polynomial g(x)

global SETTINGS

m = GetDegree(fx);
n = GetDegree(gx);


% Get Mean of entries of f(x) in C_{n-k}(f)
lambda = GetMean(fx,n-1);

% Get Mean of entries of g(x) in C_{m-k}(g)
mu = GetMean(gx,m-1);

% Divide f(x) and g(x) by geometric mean
fx_n = fx./ lambda;
gx_n = gx./ mu;

switch SETTINGS.BOOL_ALPHA_THETA
    case 'y'
        
        % Get opitmal values of alpha and theta
        [alpha, theta] = GetOptimalAlphaAndTheta(fx_n,gx_n);
        
        % Obtain f(w) and g(w) from f(x) and g(x)]
        fw = GetWithThetas(fx_n,theta);
        gw = GetWithThetas(gx_n,theta);
        
        % Plot the unprocessed and preprocessed coefficients of
        % f(x), f(w), g(x) and g(w).
        
        switch SETTINGS.PLOT_GRAPHS
            case 'y'
                PlotCoefficients(fx,fw,'f');
                PlotCoefficients(gx,gw,'g');
            case 'n'
            otherwise
                error('err')
        end
    case 'n'
        
        % Dont apply preprocessing
        theta = 1;
        alpha = 1;
        
    otherwise
        error('bool_preproc either y or n');
end

end