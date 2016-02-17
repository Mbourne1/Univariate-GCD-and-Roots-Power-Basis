function [] = o_intersection_implicit()
% given two curves f(x) and g(x), calculate their intersection

global PLOT_GRAPHS
PLOT_GRAPHS = 'y'

fx = [1; 5; 2; 1];
gx = [5; 6; 7; 8];

% Get degree of polynomial f(x)
[r,~] = size(fx);
m = r - 1;

% Get degree of polynomial g(x)
[r,~] = size(gx);
n = r - 1;

%%
% Get the polynomial h(x) = f(x) - g(x)

if m < n
    % If m > n, padd g(x) with zeros
    fx = [fx ; zeros(n-m,1)];
elseif m > n
    % If m < n, padd f(x) with zeros
    gx = [gx ; zeros(m-n,1)];
end

hx = fx-gx

%%
% Plot f(x), g(x) and h(x)

% plot f(x) and g(x)
t = linspace(-10,10,200);
f_y = polyval(fx,t);
g_y = polyval(gx,t);
h_y = polyval(hx,t);

switch PLOT_GRAPHS
    case 'y'
        figure('name','Curve Plot')
        hold on
        plot(t,f_y,'DisplayName','f(y)')
        plot(t,g_y,'DisplayName','g(y)')
        plot(t,h_y,'DisplayName','d(y)')
        legend(gca,'show');
        hold off
    case 'n'
    otherwise
        error('PLOT_GRAPHS is either y or n')
end

o_roots_mymethod(hx)
