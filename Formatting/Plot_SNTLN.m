switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s - Residuals',mfilename);
        figure('name',figure_name)
        hold on
        plot(log10(condition),'-s');
        hold off
        
    case 'n'
    otherwise
        error('err')
end