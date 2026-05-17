function [x, history, iter] = gauss_seidel_method(A, b, tol, maxit, x0)
% GAUSS_SEIDEL_METHOD Solve Ax=b by Gauss-Seidel iteration.
% Tridiagonal implementation for the 1D Poisson-type system.

x = x0;
n = length(b);

d = diag(A);
lower = diag(A, -1);
upper = diag(A, 1);

normb = norm(b);
if normb == 0
    normb = 1;
end

history = zeros(maxit + 1, 1);
history(1) = norm(b - A*x) / normb;

for iter = 1:maxit
    xold = x;

    for i = 1:n
        s = b(i);

        if i > 1
            s = s - lower(i-1) * x(i-1);
        end

        if i < n
            s = s - upper(i) * xold(i+1);
        end

        x(i) = s / d(i);
    end

    history(iter + 1) = norm(b - A*x) / normb;

    if history(iter + 1) < tol
        history = history(1:iter + 1);
        return;
    end
end
end