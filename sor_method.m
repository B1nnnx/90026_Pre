function [x, history, iter] = sor_method(A, b, omega, tol, maxit, x0)
% SOR_METHOD Solve Ax=b by successive over-relaxation.
% One-step update follows Burkardt's sor1.m.

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

        x(i) = (1 - omega) * xold(i) + omega * x(i);
    end

    history(iter + 1) = norm(b - A*x) / normb;
    if history(iter + 1) < tol
        history = history(1:iter + 1);
        return;
    end
end
end
