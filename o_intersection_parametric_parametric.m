function [] = o_intersection_parametric_parametric()
% Given two integral parametrically defined curves, obtain the points of
% intersection.

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

global BOOL_PREPROC
BOOL_PREPROC = 'y';


% fxt : vector of coefficients of x(t) with increasing powers of t.
% fyt : vector of coefficients of y(t) with increasing powers of t.
C1_xt = [0.1; 0.8; -0.1];
C1_yt = [0.5; 3; -2.5];

C2_xt = [0; 2];
C2_yt = [0; 2];

% Get the implicit representation of C_{1}(x(t),y(t))
C1_implicit = Implicitize_Integral_Parametric(C1_xt,C1_yt);

% Get implicit Representation
fprintf('The Implicit Equation of C_{1}(x,y) is given by')
fprintf('\n')
disp(C1_implicit)

PrintCoefficientsBivariate(C1_implicit,'C1')


% Get the number of rows and columns in C1(x,y)
[r,c] = size(C1_implicit);

poly = 0;

% for each row in C_{1}(x,y)
for i = 0:1:r-1
    % for each column in C_{1}(x,y)
    for j = 0:1:c-1
        
        x_component = 1;
        for k = 0:1:i
            x_component = conv(C2_xt,x_component);
        end
        
        y_component = 1;
        for k = 0:1:j
            y_component = conv(C2_yt,y_component);
        end
        
        uij = conv(x_component,y_component);
        
        uij = uij .* C1_implicit(i+1,j+1);
        
        poly = PolyAdd(poly, uij);
    end
end

fprintf('The Curve C_{3} is given by \n')
fprintf('\n')
PrintPoly(poly,'C3')

% Substitute x(t) and y(t) from C_{2} into implicit representation of C_{1}

% Get roots of C_{3} in terms of t
o_roots_mymethod(poly);
arrRoot = roots(fliplr(poly'));

fprintf('The roots in terms of t are given by')
fprintf('\n')
disp(arrRoot)

[r,~] = size(arrRoot);

x = zeros(1,r);
y = zeros(1,r);

% For each root in arrRoot
for i = 0:1:r-1
    
    % Get the ith root
    rt = arrRoot(i+1);
    
    % Substitute in to x(t)
    x_sum = 0;
    
    % for each coefficient in C1_x(t)
    for j = 0:1:size(C2_xt,1)-1
        x_sum = x_sum + (C2_xt(j+1) .* (rt.^j));
    end
    
    x(i+1) = x_sum;
    
    y_sum = 0;
    for j = 0:1:size(C2_yt) -1
        y_sum = y_sum + (C2_yt(j+1) .* (rt.^j));
    end
    
    y(i+1) = y_sum;
    
end

fprintf('The intersections are given by')
fprintf('\n')
disp([x;y])


% Substitute root values back into C_{1} x(t) and y(t) to obtain set of
% intersection points
switch PLOT_GRAPHS
    case 'y'
        
        t = linspace(-100,100,2000);
        C1_x_vals = polyval(fliplr(C1_xt'),t);
        C1_y_vals = polyval(fliplr(C1_yt'),t);
        
        C2_x_vals = polyval(fliplr(C2_xt'),t);
        C2_y_vals = polyval(fliplr(C2_yt'),t);
        
        
        C3 = polyval(fliplr(poly'),t);
        
        figure('name','Plotting C_{1} and C_{2}')
        hold on
        plot(C1_x_vals,C1_y_vals,'blue')
        plot(C2_x_vals,C2_y_vals,'red')
        plot(t,C3,'green')
        axis([-10,10,-10,10])
        grid on
        hold off
        
    case 'n'
end


end