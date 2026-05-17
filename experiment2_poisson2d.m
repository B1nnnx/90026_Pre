function results = experiment2_poisson2d()
% Compare Jacobi, Gauss-Seidel and SOR on the 2D Poisson problem.

code_dir = fileparts(mfilename('fullpath'));
addpath(code_dir);
figs_dir = fullfile(fileparts(code_dir), 'figs');
mkdir(figs_dir);

m = 20;
h = 1 / (m + 1);
N = m^2;
tol = 1e-8;
maxit = 50000;

e = ones(m, 1);
T = spdiags([-e, 2*e, -e], -1:1, m, m);
A = kron(speye(m), T) + kron(T, speye(m));
b = h^2 * ones(N, 1);
x0 = zeros(N, 1);

rho_jacobi = cos(pi / (m + 1));
omega_opt = 2 / (1 + sqrt(1 - rho_jacobi^2));

[x_j, hist_j, it_j] = jacobi_method(A, b, tol, maxit, x0);
[x_gs, hist_gs, it_gs] = gauss_seidel_method(A, b, tol, maxit, x0);
[x_sor, hist_sor, it_sor] = sor_method(A, b, omega_opt, tol, maxit, x0);

fprintf('\nExperiment 2: 2D Poisson\n');
fprintf('m = %d, N = %d\n', m, N);
fprintf('omega_opt = %.6f\n', omega_opt);
fprintf('Jacobi iterations: %d\n', it_j);
fprintf('Gauss-Seidel iterations: %d\n', it_gs);
fprintf('SOR iterations: %d\n', it_sor);

figure;
loglog(1:it_j, hist_j(2:end), 'LineWidth', 1.5); hold on;
loglog(1:it_gs, hist_gs(2:end), 'LineWidth', 1.5);
loglog(1:it_sor, hist_sor(2:end), 'LineWidth', 1.8);
grid on;
xlabel('Iteration');
ylabel('Relative residual');
title('Convergence history for 2D Poisson problem');
legend('Jacobi', 'Gauss-Seidel', sprintf('SOR, \\omega_{opt} = %.3f', omega_opt), ...
    'Location', 'southwest');
ylim([tol, 1]);
saveas(gcf, fullfile(figs_dir, 'poisson2d_convergence.pdf'));

results = struct();
results.problem = '2D Poisson';
results.m = m;
results.N = N;
results.omega_opt = omega_opt;
results.jacobi_iterations = it_j;
results.gauss_seidel_iterations = it_gs;
results.sor_iterations = it_sor;
results.x_jacobi = x_j;
results.x_gauss_seidel = x_gs;
results.x_sor = x_sor;
results.figs_dir = figs_dir;
end
