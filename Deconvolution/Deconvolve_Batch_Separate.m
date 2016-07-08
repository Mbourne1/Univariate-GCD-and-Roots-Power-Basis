function hx = Deconvolve_Batch_Separate(arr_fx)

nPolys = length(arr_fx);

for i = 1:1:nPolys-1
    hx{i} = Deconvolve(arr_fx{i},arr_fx{i+1}) ;
end

end