function fx = Examples_Roots(ex_num)


switch ex_num
    case '1'
        roots_fx = ...
            [
            0.1 1;
            0.2 2;
            0.5 3;
            ];
    case '2'   
    roots_fx = ...
            [
                6   2;
                -3  3;
                -7  1;
            ];
    case '3'
        roots_fx = ...
            [
            1   1
            2   2
            3   3
            4.5 1
            ];
end

writeToText(roots_fx,'f')

fx = get_Coeff(roots_fx);

PrintPoly(fx,'f')

end