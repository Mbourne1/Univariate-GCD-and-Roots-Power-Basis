function der_f = Differentiate(fx)
% Get the derivative of the input polynomial f(x)

% Inputs.
%
% fx :  Coefficients of polynomial f(x)
%

%%

% Get the degree of polynomial f(x)
m = GetDegree(fx);

% Get the derivative of f(x)
der_f = ((0:1:m)').*fx;

% Remove the first entry
der_f(1) = [];

end