function roots_fx = Examples_Roots(ex_num)


switch ex_num
    case '1'
        roots_fx = ...
            [
            0.1 1;
            0.2 2;
            0.5 3;
            ];
    case '1b'
        roots_fx = ...
            [
            0.1 1;
            0.2 2;
            0.9 3;
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
    case '4'
        roots_fx = ...
            [
            1.05467 1
            2.24587 2
            5.54743 3
            1.75647 2
            ];
    case '4b'
        roots_fx = ...
            [
            1.05467 1
            2.24587 2
            5.54743 3
            1.75647 2
            2.56478 5
            1.15445 1
            ];
    otherwise 
        error('err')
end




end