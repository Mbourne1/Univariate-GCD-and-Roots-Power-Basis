
switch PLOT_GRAPHS
    case 'y'
        % Plot the coefficients of fx
        figure('name','Coefficients of f(x)')
        plot(fx,'-s','DisplayName','Coefficients Before Processing')
        hold on
        plot(fw,'-o','DisplayName','Coefficients After Processing')
        legend(gca,'show');
        hold off
        
        % Plot the coefficients of g(x)
        figure('name','Coefficients of g(x)')
        plot(gx,'-s','DisplayName','Coefficients Before Processing')
        hold on
        plot(gw,'-o','DisplayName','Coefficients After Processing')
        legend(gca,'show');
        hold off
    case 'n'
    otherwise
        error('err : plot graphs either y or n')
end
