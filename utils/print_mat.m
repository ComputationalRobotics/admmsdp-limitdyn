function print_mat(A, acc)
    if nargin == 1
        acc = 0;
    end
    m = size(A, 1);
    n = size(A, 2);
    acc_print = sprintf("%%1.%de ", acc);
    for i = 1: m
        for j = 1: n
            if (A(i, j) > 0) 
                fprintf("+" + acc_print, abs(A(i, j)));
            else 
                fprintf("-" + acc_print, abs(A(i, j)));
            end
        end
        fprintf("\n");
    end
    fprintf("\n");
end

