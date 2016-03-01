function [] = o_IntersectionExplicit(ex_num_f,ex_num_g)
% given two curves f(x) and g(x), calculate their intersection

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

global BOOL_PREPROC
BOOL_PREPROC = 'n';

global BOOL_SNTLN
BOOL_SNTLN = 'n';

global MAX_ERROR_SNTLN
MAX_ERROR_SNTLN = 1e-11;

global MAX_ITE_SNTLN
MAX_ITE_SNTLN = 50;

% Get the type of example used
EXAMPLE_TYPE = 'Roots';
%EXAMPLE_TYPE = 'Coefficients';

switch EXAMPLE_TYPE
    case 'Coefficients'
        switch ex_num
            case '1'
                fx = [  972; -5076;  10503; -11502;  7415; -2916;  689; -90; 5];
                gx = [-7776; 25056; -32940;  23004; -9181;  2050; -216;   2; 1];
        end
    case 'Roots'
        % Get the coefficients of two implicitly defined polynomial curves.
        fx = Examples_Intersection_Implicit(ex_num_f);
        gx = Examples_Intersection_Implicit(ex_num_g);
        PrintPoly(fx,'f');
        PrintPoly(gx,'g');
        
        
end

% Get degree of polynomial f(x)
[r,~] = size(fx);
m = r - 1;

% Get degree of polynomial g(x)
[r,~] = size(gx);
n = r - 1;

%%

INTERSECTION_METHOD = 'Common Factors First';
%INTERSECTION_METHOD = 'All Roots';

switch INTERSECTION_METHOD
    case 'All Roots'
        % Get the polynomial h(x) = f(x) - g(x)
        % For this, f(x) and g(x) must be of the same length.
        
        % Create a zeros vector of the length of h(x)
        zero_vector = zeros(min(m,n),1);
        
        % Replace f(x) by f(x) padded with zeros for coefficients
        % x^{m+1},...
        f = zero_vector;
        f(1:m+1) = fx;
        fx_inflated = f;
        
        % Replace g(x) by g(x) padded with zeros for coefficients
        % x^{n+1},...
        g = zero_vector;
        g(1:n+1) = gx;
        gx_inflated = g;
        
        % Get the coefficients h(x)
        hx = fx_inflated - gx_inflated;
        
        % Check for zero coefficients
        if hx(end) == 0
            hx(end) = []
        end
        
        %fprintf('Roots by Matlab Method \n')
        %roots(flipud(hx))
        
    case 'Common Factors First'
        
        % Get the GCD of f(x) and g(x)
        [~,~,dx, ux, vx , alpha,theta, ~, lambda, mu] = o1(fx,gx);
        
        % Get the roots of polynomial d(x)
        my_roots = o_roots_mymethod(dx);
        
        % Get degree of u(x)
        [r,~] = size(ux);
        m = r -1;
        
        % Get degree of v(x)
        [r,~] = size(vx);
        n = r - 1;
        
        vZeros = zeros(max(m,n),1);
        u = vZeros;
        u(1:m+1) = ux;
        
        v = vZeros;
        v(1:n+1) = vx;
        
        
        % Get the coefficients h(x).
        hx = ux - vx;
        
        % Check for zero coefficients
        if hx(end) == 0
            hx(end) = [];
        end
end

% Get the degree of h(x)
[r,~] = size(hx);
m = r-1;
if m ~=0
    % Get the roots of h(x)
    roots2 = o_roots_mymethod(hx);
end


% %%
% Plot f(x), g(x) and h(x)
% plot f(x) and g(x)
t = -5:0.01:5;
f_y = polyval(flipud(fx),t);
g_y = polyval(flipud(gx),t);
h_y = polyval(flipud(hx),t);

switch PLOT_GRAPHS
    case 'y'
        figure('name','Curve Plot : f(y) and g(y)')
        hold on
        plot(t,f_y,'DisplayName','f(y)')
        plot(t,g_y,'DisplayName','g(y)')
        ylim([-1,1])
        legend(gca,'show');
        hold off
        
        figure('name','Curve Plot : f(y)-g(y)')
        hold on
        plot(t,h_y,'DisplayName','h(y)')
        ylim([-1,1])
        hold off
    case 'n'
    otherwise
        error('PLOT_GRAPHS is either y or n')
end

try
display(roots2)
catch
end


% Get the points of intersection.
% for each root in roots




end
