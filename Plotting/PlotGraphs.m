

[St,~] = dbstack();
calling_function = St(2).name;

global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Plot QR Scatter Data',mfilename)]);
        figure('name',figure_name)
        hold on
        scatter(matrix(:,1),log10(matrix(:,2)))
        hold off
        
        % Plot the minimum singular values of S_{k} k=1,...,min(m,n)
        figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Plot Min Sing Val',mfilename)]);
        figure('name',figure_name)
        hold on
        plot(log10(vMinimumSingularValues),'-o')
        plot(log10(10^-2) * ones(m+n))
        plot(log10(10^-3) * ones(m+n))
        plot(log10(10^-4) * ones(m+n))
        plot(log10(10^-5) * ones(m+n))
        plot(log10(10^-6) * ones(m+n))
        plot(log10(10^-7) * ones(m+n))
        plot(log10(10^-8) * ones(m+n))
        plot(log10(10^-9) * ones(m+n))
        plot(log10(10^-10) * ones(m+n))
        hold off
        
        % Plot the minimum distances of S_{k} k=1,...,min(m,n)
        figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Plot Min Residual',mfilename)]);
        figure('name',figure_name)
        hold on
        plot(log10(vMinimumResidual) ,'-s')
        hold off
        
        figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Max:min Row Diagonals',mfilename)]);
        figure('name',figure_name)
        
        vRatio_MaxMin_Diagonals_R = vMaxDiagR1 ./ vMinDiagR1;
        plot(log10(vRatio_MaxMin_Diagonals_R),'red-s');
        hold on
        axis([1,min(m,n),0,inf])
        legend('Max:Min diag element of subresultant S_{k}');
        title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
        ylabel('log_{10} max:min diag element')
        hold off
        
        % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
        figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Max:min Row Norms',mfilename)]);
        figure('name',figure_name)
        
        vRatio_MaxMin_RowNorm_R = vMaxRowNormR1./vMinRowNormR1;
        plot(log10(vRatio_MaxMin_RowNorm_R),'red-s');
        hold on
        axis([1,min(m,n),0,inf])
        legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
        title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
        hold off
        
        
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
        
end