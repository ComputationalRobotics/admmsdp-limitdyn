function y = chol_solve(rhsy, R, P)
    rhsy = P' * rhsy;
    tmp = R' \ rhsy;
    tmp = R \ tmp;
    y = P * tmp;
end