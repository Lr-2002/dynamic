clc,clear;

addpath("dyn_iden_fun");

urdf_model_dir =  "/home/zhangdeyong/ros2_ws/tora_one_dyn_drag/toraone_dyn_iden/pr_tora_one_dev_description/urdf/pr_tora_one_dev.urdf";

robot = importrobot(urdf_model_dir);

s=randomConfiguration(robot);

% 获取各个刚体
L{1} = robot.getBody('larm_1_Link');
L{2} = robot.getBody('larm_2_Link');
L{3} = robot.getBody('larm_3_Link');
L{4} = robot.getBody('larm_4_Link');
L{5} = robot.getBody('larm_5_Link');
L{6} = robot.getBody('larm_6_Link');
L{7} = robot.getBody('larm_7_Link');

% 读取质量
M1 = L{1}.Mass;                  
M2 = L{2}.Mass; 
M3 = L{3}.Mass; 
M4 = L{4}.Mass; 
M5 = L{5}.Mass;
M6 = L{6}.Mass;
M7 = L{7}.Mass;
M = [M1, M2, M3, M4, M5, M6, M7];

% 读取重心位置
Pc{1} = L{1}.CenterOfMass';
Pc{2} = L{2}.CenterOfMass';
Pc{3} = L{3}.CenterOfMass'; 
Pc{4} = L{4}.CenterOfMass'; 
Pc{5} = L{5}.CenterOfMass';
Pc{6} = L{6}.CenterOfMass'; 
Pc{7} = L{7}.CenterOfMass'; 

% 读取惯量 [Ixx, Iyy, Izz, Iyz, Ixz, Ixy]
Ic{1} = getInertiaMatrix( L{1}.Inertia );
Ic{2} = getInertiaMatrix( L{2}.Inertia );
Ic{3} = getInertiaMatrix( L{3}.Inertia );
Ic{4} = getInertiaMatrix( L{4}.Inertia );
Ic{5} = getInertiaMatrix( L{5}.Inertia );
Ic{6} = getInertiaMatrix( L{6}.Inertia );
Ic{7} = getInertiaMatrix( L{7}.Inertia );

% 关节位置，速度，加速度
% q = [0, -0.5, 0, -1.57, 0, 0, 0];
q = [0, 0, 0, 0, 0, 0, 0];
dq = [0, 0, 0, 0, 0, 0, 0];
ddq = [0, 0, 0, 0, 0, 0, 0];

% 牛顿欧拉法计算
tau = tau_newtonEuler(q,dq,ddq,M,Ic,Pc)


robot.show(s)

% 导入urdf模型到simuscape，生成simulink模型文件
% smimport(urdf_model_dir)

