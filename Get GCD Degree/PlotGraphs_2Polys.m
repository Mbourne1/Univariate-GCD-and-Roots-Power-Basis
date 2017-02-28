

[St,~] = dbstack();
calling_function = St(2).name;

global SETTINGS
if( SETTINGS.PLOT_GRAPHS)
    
        
        plot_qr = false;
        plot_minimumSingularValues = true;
        plot_minDistance = false;
        plot_MaxMinDiags = false;
        plot_MaxMinRowNorm = false;
        
        if plot_qr == true
            figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Plot QR Scatter Data',mfilename)]);
            figure('name',figure_name)
            xlim([1 min(m,n)]);
            hold on
            scatter(matrix(:,1),log10(matrix(:,2)))
            vline([lower_lim upper_lim],{'b','b'},{'',''});
            hold off
        end
        
        if plot_minimumSingularValues == true
            % Plot the minimum singular values of S_{k} k=1,...,min(m,n)
            figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Plot Min Sing Val',mfilename)]);
            figure('name',figure_name)
            hold on
            xlim([1 min(m,n)]);
            plot(log10(vMinimumSingularValues),'-o')
            vline([lower_lim upper_lim],{'b','b'},{'',''});
            hold off
        end
        
        if plot_minDistance == true
            % Plot the minimum distances of S_{k} k=1,...,min(m,n)
            figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Plot Min Residual',mfilename)]);
            figure('name',figure_name)
            hold on
            xlim([1 min(m,n)]);
            plot(log10(vMinimumResidual) ,'-s')
            vline([lower_lim upper_lim],{'b','b'},{'',''});
            hold off
        end
        
        
        if plot_MaxMinDiags == true
            figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Max:min Row Diagonals',mfilename)]);
            figure('name',figure_name)
            vRatio_MaxMin_Diagonals_R = vMaxDiagR1 ./ vMinDiagR1;
            plot(log10(vRatio_MaxMin_Diagonals_R),'red-s');
            hold on
            xlim([1 min(m,n)]);
            legend('Max:Min diag element of subresultant S_{k}');
            title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
            ylabel('log_{10} max:min diag element')
            vline([lower_lim upper_lim],{'b','b'},{'',''});
            hold off
        end
        
        if plot_MaxMinRowNorm == true
            % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
            figure_name = sprintf([calling_function ' : ' mfilename ' : ' sprintf('%s : Max:min Row Norms',mfilename)]);
            figure('name',figure_name)

            vRatio_MaxMin_RowNorm_R = vMaxRowNormR1./vMinRowNormR1;
            plot(log10(vRatio_MaxMin_RowNorm_R),'red-s');
            hold on
            xlim([1 min(m,n)]);
            legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
            title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
            vline([lower_lim upper_lim],{'b','b'},{'',''});
            hold off
        end
        
    
        
end