# 本matlab仓库主要用于动力学参数辨识中，激励轨迹优化、生成、最小参数集辨识计算

# 1. 激励轨迹优化matlab脚本
运行TrajPSOID_8dof.m，会在工作空间生成x_opt变量，保存成x_opt.mat

# 1b. 6-DOF 激励轨迹优化matlab脚本
运行TrajPSOID_6dof.m，会在工作空间生成x_opt变量，保存成x_opt_6dof.mat

# 2. 激励轨迹生成matlab脚本
运行pot_final_trajectory_script.m，会将轨迹保存在traj.txt文件里

# 2b. 6-DOF 激励轨迹生成matlab脚本
运行pot_final_trajectory_script_6dof.m，会将轨迹保存在traj_6dof.txt文件里（末尾补零到8轴）

# 3. 动力学参数辨识matlab脚本
运行Pr_iden_arm_real.m，P_min_dh即为最小参数集

# 3b. 6-DOF 动力学参数辨识matlab脚本
运行Pr_iden_arm_real_6dof.m，P_min_dh_6dof即为6轴最小参数集
