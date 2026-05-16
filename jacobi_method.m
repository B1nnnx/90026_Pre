function [x, history, iter] = jacobi_method(A, b, tol, maxit, x0)
% JACOBI_METHOD Solve Ax=b by Jacobi iteration.
% One-step update follows Burkardt's jacobi3.m:
% x_new = x_old + D^{-1}(b - A*x_old).

x = x0;
D = diag(A);
normb = norm(b);

history = zeros(maxit + 1, 1);
history(1) = norm(b - A*x) / normb;

for iter = 1:maxit
    xold = x;
    x = xold + (b - A*xold) ./ D;

    history(iter + 1) = norm(b - A*x) / normb;
    if history(iter + 1) < tol
        history = history(1:iter + 1);
        return;
    end
end
end
