function [M, rho, xi] = util_ocsvm(gam, nu)
m = size(gam, 1); gsize = size(gam, 2) - 1; a = 1/(nu*m); 
cvx_begin
    variable rho nonnegative;
    variable M(gsize, gsize) semidefinite;
    variable xi(m) nonnegative;
    minimize(sum(sum_square(M))/2 - rho + a*sum(xi));
    for i = 1:m
        trace(M * reshape(gam(i, 2:end, 2:end), gsize, gsize)) - 1*gam(i, 1, 1) >= rho - xi(i); 
    end
cvx_end
end