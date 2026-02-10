clc, clear

% Trajectory settings
n_joints = 6;
T = 30;
dt = 0.002;
t = 0:dt:T;
N = 5;
freq_base = 1/T;

x_opt_file = 'x_opt_6dof.mat';
load(x_opt_file);

[q, dq, ddq] = plot_final_trajectory(x_opt, n_joints, N, t, freq_base);

% Pad to 8 joints for downstream tools
q_out = [q; zeros(2, size(q,2))];
out_file = 'traj_6dof.txt';
writematrix(q_out', out_file, 'Delimiter', ' ');

%% Final trajectory plotting
function [q, dq, ddq] = plot_final_trajectory(x_opt, n_joints, N, t, freq_base)
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
