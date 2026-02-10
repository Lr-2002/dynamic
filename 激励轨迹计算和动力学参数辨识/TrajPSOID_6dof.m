%% Fourier excitation trajectory - PSO (6-DOF)
clc; clear; close all;

addpath('dyn_iden_fun/');

% ====================== 1. Globals ======================
global iteration_history best_cost_history best_solution_history traj_plots;
global n_joints N t freq_base;
iteration_history = [];
best_cost_history = [];
best_solution_history = [];
traj_plots = cell(6,3);

% Trajectory settings
n_joints = 6;
T = 30;
dt = 0.05;
t = 0:dt:T;
N = 5;
freq_base = 1/T;

% Joint limits (from URDF, first 6 joints)
q_min = [-0.72, -1.29, -1.81, -1.538, -0.48, -1.481];
q_max = [ 0.714,  0.5,   0.56,  1.5329, 1.57,  1.458];

% Velocity limits (rad/s)
dq_max = [70, 40, 40, 36, 36, 120] * 2*pi / 60 * 0.7;
ddq_max = dq_max * 2;

% PSO options
pso_options = optimoptions('particleswarm', ...
    'SwarmSize', 30, ...
    'MaxIterations', 10000, ...
    'Display', 'iter', ...
    'UseVectorized', false, ...
    'FunctionTolerance', 1e-6, ...
    'SelfAdjustmentWeight', 1.49, ...
    'SocialAdjustmentWeight', 1.49, ...
    'OutputFcn', @pso_outfun);

% Variable bounds
n_vars = n_joints * (2*N + 1);
lb = -3 * ones(1, n_vars);
ub = 3 * ones(1, n_vars);

% ====================== 2. Visualization ======================
figure('Position', [100, 100, 1400, 800]);
subplot(3,4,1); hold on; grid on;
title('Convergence');
xlabel('Iteration'); ylabel('Cost');
set(gca, 'YScale', 'log');
convergence_plot = plot(NaN, NaN, 'b-', 'LineWidth', 2);

for i = 1:n_joints
    subplot(3,4,i+1); hold on; grid on;
    traj_plots{i,1} = plot(t, NaN(size(t)), 'b-', 'LineWidth', 2);
    traj_plots{i,2} = plot(t, NaN(size(t)), 'r--', 'LineWidth', 1.5);
    traj_plots{i,3} = plot(t, NaN(size(t)), 'g-.', 'LineWidth', 1.5);
    title(sprintf('Joint %d (iter: 0)', i));
    xlabel('Time (s)'); ylabel('Value');
    legend('pos', 'vel', 'acc', 'Location', 'best');
    ylim([-5, 5]);
end

% ====================== 3. PSO ======================
fprintf('Starting PSO...\n');
fun = @(x) MYCOND_PSO(x, n_joints, N, q_min, q_max, dq_max, ddq_max, t, freq_base);
[x_opt, fval, exitflag, output] = particleswarm(fun, n_vars, lb, ub, pso_options);

% ====================== 4. Results ======================
fprintf('Done.\n');
disp(['Best cost: ', num2str(fval)]);
disp(['Iterations: ', num2str(output.iterations)]);

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1);
plot(iteration_history, best_cost_history, 'b-', 'LineWidth', 2);
grid on; title('Convergence'); xlabel('Iteration'); ylabel('Cost');
set(gca, 'YScale', 'log');

subplot(1,2,2);
semilogy(iteration_history, best_cost_history, 'r-', 'LineWidth', 2);
grid on; title('Convergence (log)'); xlabel('Iteration'); ylabel('Cost');

plot_final_trajectory(x_opt, n_joints, N, t, freq_base);

%% ====================== Subfunctions ======================
function stop = pso_outfun(optimValues, state)
    global iteration_history best_cost_history best_solution_history n_joints N t freq_base;

    stop = false;

    switch state
        case 'init'
            fprintf('PSO init done\n');
        case 'iter'
            iteration = optimValues.iteration;
            best_cost = optimValues.bestfval;
            best_solution = optimValues.bestx;

            iteration_history(end+1) = iteration;
            best_cost_history(end+1) = best_cost;
            best_solution_history(end+1, :) = best_solution;

            subplot(3,4,1);
            plot(iteration_history, best_cost_history, 'b-', 'LineWidth', 2);
            title(sprintf('Convergence (iter: %d, best: %.4f)', iteration, best_cost));
            drawnow;

            if mod(iteration, 5) == 0 || iteration == 1
                update_trajectory_plots(best_solution, iteration, n_joints, N, t, freq_base);
            end

            fprintf('Iter %d: best = %.6f\n', iteration, best_cost);
        case 'done'
            fprintf('PSO finished\n');
    end
end

function update_trajectory_plots(x, iteration, n_joints, N, t, freq_base)
    global traj_plots;

    a = zeros(n_joints, N);
    b = zeros(n_joints, N);
    qs = zeros(n_joints, 1);

    for i = 1:n_joints
        for j = 1:N
            a(i, j) = x((i-1)*(2*N+1) + j);
            b(i, j) = x((i-1)*(2*N+1) + N + j);
        end
        qs(i) = x((i-1)*(2*N+1) + 2*N + 1);
    end

    q = zeros(n_joints, length(t));
    dq = zeros(n_joints, length(t));
    ddq = zeros(n_joints, length(t));

    for j = 1:length(t)
        for i = 1:n_joints
            sum_q = 0; sum_dq = 0; sum_ddq = 0;
            for k = 1:N
                w_k = 2*pi * k * freq_base;
                sum_q = sum_q + a(i,k)/w_k * sin(w_k*t(j)) - b(i,k)/w_k * cos(w_k*t(j));
                sum_dq = sum_dq + a(i,k) * cos(w_k*t(j)) + b(i,k) * sin(w_k*t(j));
                sum_ddq = sum_ddq -a(i,k)*w_k * sin(w_k*t(j)) + b(i,k)*w_k * cos(w_k*t(j));
            end
            q(i,j) = qs(i) + sum_q;
            dq(i,j) = sum_dq;
            ddq(i,j) = sum_ddq;
        end
    end
    q = q - q(:, 1);

    for i = 1:n_joints
        subplot(3,4,i+1);
        set(traj_plots{i,1}, 'XData', t, 'YData', q(i,:));
        set(traj_plots{i,2}, 'XData', t, 'YData', dq(i,:));
        set(traj_plots{i,3}, 'XData', t, 'YData', ddq(i,:));
        title(sprintf('Joint %d (iter: %d)', i, iteration));
    end

    drawnow;
end

function cost = MYCOND_PSO(x, n_joints, N, q_min, q_max, dq_max, ddq_max, t, freq_base)

    a = zeros(n_joints, N);
    b = zeros(n_joints, N);
    qs = zeros(n_joints, 1);

    for i = 1:n_joints
        for j = 1:N
            a(i, j) = x((i-1)*(2*N+1) + j);
            b(i, j) = x((i-1)*(2*N+1) + N + j);
        end
        qs(i) = x((i-1)*(2*N+1) + 2*N + 1);
    end

    [base_cost, cond_cost] = calculate_condition(a, b, qs);
    penalty = calculate_constraint_penalty(a, b, qs, q_min, q_max, dq_max, ddq_max, t, freq_base);

    cost = cond_cost + 1e3*penalty;
end

function penalty = calculate_constraint_penalty(a, b, qs, q_min, q_max, dq_max, ddq_max, t, freq_base)
    penalty = 0;
    n_joints = size(a,1);
    N = size(a,2);

    penalty = penalty + sum(max(0, abs(a(:)) - 3).^2);
    penalty = penalty + sum(max(0, abs(b(:)) - 3).^2);

    [q, dq, ddq] = calculate_trajectory(a, b, qs, t, freq_base, N, n_joints);

    for i = 1:n_joints
        penalty = penalty + sum(max(0, q(i,:)-q_max(i)) + max(0, q_min(i)-q(i,:)));
    end
    for i = 1:n_joints
        penalty = penalty + sum(max(0, abs(dq(i,:))-dq_max(i)));
    end
    for i = 1:n_joints
        penalty = penalty + sum(max(0, abs(ddq(i,:))-ddq_max(i)));
    end

    penalty = penalty + sum(abs(dq(:,1)));
    penalty = penalty + sum(abs(ddq(:,1)));
end

function plot_final_trajectory(x_opt, n_joints, N, t, freq_base)
    figure('Position', [100, 100, 1200, 600]);

    a = zeros(n_joints, N);
    b = zeros(n_joints, N);
    qs = zeros(n_joints, 1);

    for i = 1:n_joints
        for j = 1:N
            a(i, j) = x_opt((i-1)*(2*N+1) + j);
            b(i, j) = x_opt((i-1)*(2*N+1) + N + j);
        end
        qs(i) = x_opt((i-1)*(2*N+1) + 2*N + 1);
    end

    [q, dq, ddq] = calculate_trajectory(a, b, qs, t, freq_base, N, n_joints);

    for i = 1:n_joints
        subplot(2,3,i);
        plot(t, q(i,:), 'b-', 'LineWidth', 2); hold on;
        plot(t, dq(i,:), 'r--', 'LineWidth', 1.5);
        plot(t, ddq(i,:), 'g-.', 'LineWidth', 1.5);
        title(sprintf('Joint %d', i));
        xlabel('Time (s)'); ylabel('Value');
        legend('pos', 'vel', 'acc');
        grid on;
    end
end

function [q, dq, ddq] = calculate_trajectory(a, b, qs, t, freq_base, N, n_joints)
    q = zeros(n_joints, length(t));
    dq = zeros(n_joints, length(t));
    ddq = zeros(n_joints, length(t));

    for j = 1:length(t)
        for i = 1:n_joints
            sum_q = 0; sum_dq = 0; sum_ddq = 0;
            for k = 1:N
                w_k = 2*pi * k * freq_base;
                sum_q = sum_q + a(i,k)/w_k * sin(w_k*t(j)) - b(i,k)/w_k * cos(w_k*t(j));
                sum_dq = sum_dq + a(i,k) * cos(w_k*t(j)) + b(i,k) * sin(w_k*t(j));
                sum_ddq = sum_ddq -a(i,k)*w_k * sin(w_k*t(j)) + b(i,k)*w_k * cos(w_k*t(j));
            end
            q(i,j) = qs(i) + sum_q;
            dq(i,j) = sum_dq;
            ddq(i,j) = sum_ddq;
        end
    end
    q = q - q(:, 1);
end

function [cost, cond_cost] = calculate_condition(a, b, qs)
    cost = 0.1 * (sum(a(:).^2) + sum(b(:).^2)) + 0.01 * sum(qs.^2);
    x = [a, b, qs];
    cond_cost = MYCOND_6dof(x);
end
