clc,clear

% 轨迹参数
n_joints = 8;           % 关节数
T = 30;                 % 周期时间 (s)
dt = 0.002;             % 采样时间间隔
t = 0:dt:T;             % 时间序列
N = 5;                  % 傅里叶级数阶数
freq_base = 1/T;        % 基频 (Hz)

load('x_opt_2.mat');


% 显示最佳解的轨迹
[q, dq, ddq] = plot_final_trajectory(x_opt, n_joints, N, t, freq_base);

% 保存数据到txt文件
writematrix(q', 'traj.txt', 'Delimiter', ' ');


%% 最终轨迹绘制
function [q, dq, ddq] = plot_final_trajectory(x_opt, n_joints, N, t, freq_base)
    figure('Position', [100, 100, 1200, 600]);
    
    % 解析最优解
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
    
    % 计算最终轨迹
    [q, dq, ddq] = calculate_trajectory(a, b, qs, t, freq_base, N, n_joints);
    
    % 绘制最终轨迹（3×2布局，适配6个关节）
    for i = 1:n_joints
        subplot(3,3,i);
        plot(t, q(i,:), 'b-', 'LineWidth', 2); hold on;
        plot(t, dq(i,:), 'r--', 'LineWidth', 1.5);
        plot(t, ddq(i,:), 'g-.', 'LineWidth', 1.5);
        title(sprintf('关节 %d - 最终轨迹', i));
        xlabel('时间 (s)'); ylabel('幅值');
        legend('位置', '速度', '加速度');
        grid on;
    end
end

%% 轨迹计算函数
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
