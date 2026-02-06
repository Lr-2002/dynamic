%% 傅里叶激励轨迹生成 - 粒子群优化版本（带迭代过程可视化）
clc; clear; close all;


addpath('dyn_iden_fun/');

% ====================== 1. 初始化全局变量 ======================
global iteration_history best_cost_history best_solution_history traj_plots;
global n_joints N t freq_base;  % 新增：将轨迹参数声明为全局变量
iteration_history = [];       % 初始化：空数值数组（迭代次数）
best_cost_history = [];       % 初始化：空数值数组（最优目标函数值）
best_solution_history = [];   % 初始化：空数组（最优解）
traj_plots = cell(8,3);       % 初始化：8关节×3曲线（位置/速度/加速度）

% 轨迹参数
n_joints = 8;           % 关节数
T = 30;                 % 周期时间 (s)
dt = 0.05;             % 采样时间间隔
t = 0:dt:T;             % 时间序列
N = 5;                  % 傅里叶级数阶数
freq_base = 1/T;        % 基频 (Hz)

% 关节运动范围限制（从 dataarm.urdf 硬编码）
% j1..j8 limits from /home/lr-2002/project/DataArm/dataarm/urdfs/dataarm/urdf/dataarm.urdf
q_min = [-0.72, -1.29, -1.81, -1.538, -0.48, -1.481, -1.2, -1.2];
q_max = [ 0.714,  0.5,   0.56,  1.5329, 1.57,  1.458,  1.2,  1.2];

% URDF 速度上限为0（未填）。这里保留默认经验值，避免轨迹被全部惩罚为不可行
dq_max = [70, 40, 40, 36, 36, 120, 45, 45] * 2*pi / 60 * 0.7;   % 最大速度限制 (rad/s)
ddq_max = dq_max * 2;  % 最大加速度限制 (rad/s^2)

% PSO参数设置
pso_options = optimoptions('particleswarm', ...
    'SwarmSize', 30, ...          % 粒子数量
    'MaxIterations', 10000, ...      % 最大迭代次数
    'Display', 'iter', ...        % 显示迭代过程
    'UseVectorized', false, ...   % 不使用向量化
    'FunctionTolerance', 1e-6, ...% 函数值容忍度
    'SelfAdjustmentWeight', 1.49, ...
    'SocialAdjustmentWeight', 1.49, ...
    'OutputFcn', @pso_outfun);    % 添加输出函数用于记录迭代过程

% 定义变量范围和初始值
n_vars = n_joints * (2*N + 1);  % 总变量数
lb = -3 * ones(1, n_vars);      % 下界
ub = 3 * ones(1, n_vars);       % 上界

% ====================== 2. 初始化可视化界面（修复子图布局！）======================
figure('Position', [100, 100, 1400, 800]);
% 子图1：收敛曲线（3行4列布局，第1个子图）
subplot(3,4,1); hold on; grid on;
title('优化收敛曲线');
xlabel('迭代次数');
ylabel('目标函数值');
set(gca, 'YScale', 'log');
convergence_plot = plot(NaN, NaN, 'b-', 'LineWidth', 2);

% 子图2-7：6个关节的位置/速度/加速度轨迹（3行4列，索引2~7）
for i = 1:n_joints
    subplot(3,4,i+1); hold on; grid on;  % 索引2~7，适配3×4布局
    % 初始化三条曲线：位置(蓝)、速度(红)、加速度(绿)
    traj_plots{i,1} = plot(t, NaN(size(t)), 'b-', 'LineWidth', 2);   % 位置
    traj_plots{i,2} = plot(t, NaN(size(t)), 'r--', 'LineWidth', 1.5); % 速度
    traj_plots{i,3} = plot(t, NaN(size(t)), 'g-.', 'LineWidth', 1.5); % 加速度
    title(sprintf('关节 %d 轨迹 (迭代: 0)', i));
    xlabel('时间 (s)'); ylabel('幅值');
    legend('位置', '速度', '加速度', 'Location', 'best');
    ylim([-5, 5]); % 固定显示范围
end

% ====================== 3. 启动PSO优化 ======================
fprintf('开始粒子群优化...\n');
fun = @(x) MYCOND_PSO(x, n_joints, N, q_min, q_max, dq_max, ddq_max, t, freq_base);
[x_opt, fval, exitflag, output] = particleswarm(fun, n_vars, lb, ub, pso_options);

% ====================== 4. 最终结果展示 ======================
fprintf('优化完成!\n');
disp(['最优目标函数值: ', num2str(fval)]);
disp(['迭代次数: ', num2str(output.iterations)]);

% 绘制最终收敛曲线
figure('Position', [100, 100, 800, 400]);
subplot(1,2,1);
plot(iteration_history, best_cost_history, 'b-', 'LineWidth', 2);
grid on; title('优化收敛曲线'); xlabel('迭代次数'); ylabel('目标函数值');
set(gca, 'YScale', 'log');

subplot(1,2,2);
semilogy(iteration_history, best_cost_history, 'r-', 'LineWidth', 2);
grid on; title('对数尺度收敛曲线'); xlabel('迭代次数'); ylabel('目标函数值(log)');

% 显示最佳解的轨迹
plot_final_trajectory(x_opt, n_joints, N, t, freq_base);

%% ====================== 子函数定义 ======================
%% PSO输出函数 - 记录迭代过程并更新可视化
function stop = pso_outfun(optimValues, state)
    global iteration_history best_cost_history best_solution_history n_joints N t freq_base;
    
    stop = false;
    
    switch state
        case 'init'
            fprintf('PSO初始化完成\n');
            
        case 'iter'
            % 记录当前迭代信息（确保是数值类型）
            iteration = optimValues.iteration;          % 迭代次数（数值）
            best_cost = optimValues.bestfval;           % 最优目标函数值（数值）
            best_solution = optimValues.bestx;          % 最优解（数组）
            
            % 赋值（此时变量已初始化，类型匹配）
            iteration_history(end+1) = iteration;
            best_cost_history(end+1) = best_cost;
            best_solution_history(end+1, :) = best_solution;
            
            % 更新收敛曲线（对应3×4布局的第1个子图）
            subplot(3,4,1);
            plot(iteration_history, best_cost_history, 'b-', 'LineWidth', 2);
            title(sprintf('优化收敛曲线 (迭代: %d, 最优值: %.4f)', iteration, best_cost));
            drawnow; % 强制刷新图像
            
            % 每5次迭代更新一次轨迹显示
            if mod(iteration, 5) == 0 || iteration == 1
                update_trajectory_plots(best_solution, iteration, n_joints, N, t, freq_base);
            end
            
            fprintf('迭代 %d: 最优值 = %.6f\n', iteration, best_cost);
            
        case 'done'
            fprintf('PSO优化完成\n');
    end
end

%% 更新轨迹图（修复参数传递+traj_plots调用）
function update_trajectory_plots(x, iteration, n_joints, N, t, freq_base)
    global traj_plots;
    
    % 解析当前最优解
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
    
    
    % 计算轨迹
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
    
    % 更新图形（对应3×4布局的索引2~7）
    for i = 1:n_joints
        subplot(3,4,i+1);
        set(traj_plots{i,1}, 'XData', t, 'YData', q(i,:));   % 位置
        set(traj_plots{i,2}, 'XData', t, 'YData', dq(i,:));  % 速度
        set(traj_plots{i,3}, 'XData', t, 'YData', ddq(i,:)); % 加速度
        title(sprintf('关节 %d 轨迹 (迭代: %d)', i, iteration));
    end
    
    drawnow; % 强制刷新
end

%% 目标函数（包含约束惩罚，修复参数传递）
function cost = MYCOND_PSO(x, n_joints, N, q_min, q_max, dq_max, ddq_max, t, freq_base)

    % 解析变量
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
    
    % 计算基本成本
    [base_cost, cond_cost] = calculate_condition(a, b, qs);
    
    % 计算约束违反惩罚
    penalty = calculate_constraint_penalty(a, b, qs, q_min, q_max, dq_max, ddq_max, t, freq_base);
    
    % 总成本 = 基本成本 + 惩罚项
%     cost = base_cost*1e5 + cond_cost + 1e4*penalty; % 惩罚系数
    cost = cond_cost + 1e3*penalty; % 惩罚系数

%     sprintf('base_cost %f , cond_cost %f,penalty: %f', base_cost, cond_cost, penalty)
end

%% 约束惩罚计算（补充轨迹约束检查）
function penalty = calculate_constraint_penalty(a, b, qs, q_min, q_max, dq_max, ddq_max, t, freq_base)
    penalty = 0;
    n_joints = size(a,1);
    N = size(a,2);
    
    % 1. 检查系数范围
    penalty = penalty + sum(max(0, abs(a(:)) - 3).^2);
    penalty = penalty + sum(max(0, abs(b(:)) - 3).^2);
%     penalty = penalty + sum(max(0, abs(qs) - 2).^2);
    
    % 2. 计算轨迹并检查位置/速度约束（补充实际约束）
    [q, dq, ddq] = calculate_trajectory(a, b, qs, t, freq_base, N, n_joints);
    % 位置越界惩罚
    for i = 1:n_joints
        penalty = penalty + sum(max(0, q(i,:)-q_max(i)) + max(0, q_min(i)-q(i,:)));
    end
    % 速度越界惩罚
    for i = 1:n_joints
        penalty = penalty + sum(max(0, abs(dq(i,:))-dq_max(i)));
    end
    % 加速度越界惩罚
    for i = 1:n_joints
        penalty = penalty + sum(max(0, abs(ddq(i,:))-ddq_max(i)));
    end

    % 速度,加速度初始值不为0惩罚
    penalty = penalty + sum(abs(dq(:,1)));
    penalty = penalty + sum(abs(ddq(:,1)));

end

%% 最终轨迹绘制
function plot_final_trajectory(x_opt, n_joints, N, t, freq_base)
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

%% 目标函数计算（示例）
function [cost, cond_cost] = calculate_condition(a, b, qs)
    % 示例目标函数：最小化运动幅值同时保持多样性
    cost = 0.1 * (sum(a(:).^2) + sum(b(:).^2)) + 0.01 * sum(qs.^2);
    x = [a, b, qs];
    cond_cost=MYCOND_8dof(x);

end
