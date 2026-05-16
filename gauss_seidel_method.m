function [x, history, iter] = gauss_seidel_method(A, b, tol, maxit, x0)
% GAUSS_SEIDEL_METHOD Solve Ax=b by Gauss-Seidel iteration.
% One-step update follows Burkardt's gauss_seidel1.m.

x = x0;
n = length(b);
normb = norm(b);

history = zeros(maxit + 1, 1);
history(1) = norm(b - A*x) / normb;

for iter = 1:maxit
    xold = x;

    for i = 1:n
        x(i) = b(i);
        x(i) = x(i) - A(i, 1:i-1) * x(1:i-1);
        x(i) = x(i) - A(i, i+1:n) * xold(i+1:n);
        x(i) = x(i) / A(i, i);
    end

    history(iter + 1) = norm(b - A*x) / normb;
    if history(iter + 1) < tol
        history = history(1:iter + 1);
        return;
    end
end
end
