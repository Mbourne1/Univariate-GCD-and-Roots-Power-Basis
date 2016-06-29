global PLOT_GRAPHS

%% Plot data
switch PLOT_GRAPHS
    case 'y'
        
        x = lower_lim:1:upper_lim;
        
        figure_name = sprintf('%s : Plot QR Scatter Data',mfilename);
        figure('name',figure_name)
        hold on
        scatter(matrix(:,1),log10(matrix(:,2)))
        hold off
        
        % Plot the minimum singular values of S_{k} k=1,...,min(m,n)
        figure_name = sprintf('%s : Plot Min Sing Val',mfilename);
        figure('name',figure_name)
        hold on
        plot(x,log10(vMinimumSingularValues),'-o')
        hold off
        
        % Plot the minimum distances of S_{k} k=1,...,min(m,n)
        figure_name = sprintf('%s : Plot Min Distance',mfilename);
        figure('name',figure_name)
        hold on
        plot(x,log10(vMinimumDistance) ,'-s')
        hold off
        
        figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
        figure('name',figure_name)
        
        plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
        hold on
        axis([lower_lim,upper_lim,0,inf])
        legend('Max:Min diag element of subresultant S_{k}');
        title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
        ylabel('log_{10} max:min diag element')
        hold off
        
        % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Norms',mfilename);
        figure('name',figure_name)
        
        plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
        hold on
        axis([lower_lim,upper_lim,0,inf])
        legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
        title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
        hold off
        
        
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
        
end