function arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta)
% Given an array of polynomials h(\omega), get 

% Remove thetas from h_{i}(w) to get h_{i}(x)
arr_hx = cell(nPolys_arr_hx,1);
for i = 1:1:nPolys_arr_hx
    arr_hx{i} = GetWithoutThetas(arr_hw{i}, theta);
end

end

