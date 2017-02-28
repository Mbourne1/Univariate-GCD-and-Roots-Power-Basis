function [] = o_gcd_batch_Unvariate_2Polys()
% Perform a batch of gcd examples
%
%
% % Example
%
% o_gcd_batch() 

ex_num_arr = {'1','2','3','4','5'};
emin_arr = {1e-8,1e-10,1e-12};
%emin_arr = {1e-8};
emax_arr = {1e-12};
mean_method_arr = {'Geometric Mean Matlab Method'};
bool_alpha_theta_arr = {'y','n'};
low_rank_approx_method_arr = {'Standard STLN','Standard SNTLN','None'};
%apf_method_arr = {'None','Standard APF Nonlinear','Standard APF Linear'};
apf_method_arr = {'None'};


parfor i1 = 1:1:length(ex_num_arr)
    
    ex_num = ex_num_arr{i1};
    
    for i2 = 1:1:length(emin_arr)
        
        emin = emin_arr{i2};
        
        for i3 = 1:1:length(emax_arr)
            
            emax = emax_arr{i3};
            
            for i4 = 1:1:length(mean_method_arr)
                
                mean_method = mean_method_arr{i4};
                
                for i5 = 1:1:length(bool_alpha_theta_arr)
                    
                    bool_alpha_theta = bool_alpha_theta_arr{i5};
                    
                    for i6 = 1:1:length(low_rank_approx_method_arr)
                        
                        low_rank_approx_method = low_rank_approx_method_arr{i6};
                        
                        for i7 = 1:1:length(apf_method_arr)
                            
                            apf_method = apf_method_arr{i7};
                            try
                                close all;
                                clc;
                                o_gcd_2Polys(ex_num,emin,emax,mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
                                fileId = fopen('log.txt','a')
                                fprintf(fileId,'%s','success \n');
                                fclose(fileId);
                            catch err
                                fileId = fopen('log.txt','a')
                                fprintf(fileId,'%s \n\n\n',getReport(err));
                                fclose(fileId);
                            end
                            
                        end
                    end
                end
            end
        end
    end
end
