            
addpath("dyn_iden_fun/");

% 导入机械臂实际采集的关节位置，角速度，电流数据
data = importdata("traj_1_modify/forward_modify/traj_12-19-18-09-47/uts_traj_left.txt");

Y = [];
tau = [];
P_min_tar = [];

si = 1;
ei = size(data,1);
step = 1;
numOfData = (ei-si)/step+1;
dlta_t_msrd = 0.002;

tau_real = [];
tau_model = [];
tau_model_ = [];


[qf, dqf, ddqf, tf, q_rm, dq_rm, t_rm] = data_parse(data,dlta_t_msrd);
qf1 = qf{1};  qf2 = qf{2}; qf3 = qf{3}; qf4 = qf{4};
qf5 = qf{5};  qf6 = qf{6}; qf7 = qf{7}; qf8 = qf{8};
dqf1 = dqf{1};  dqf2 = dqf{2}; dqf3 = dqf{3}; dqf4 = dqf{4};
dqf5 = dqf{5};  dqf6 = dqf{6}; dqf7 = dqf{7}; dqf8 = dqf{8};
ddqf1 = ddqf{1};  ddqf2 = ddqf{2}; ddqf3 = ddqf{3}; ddqf4 = ddqf{4};
ddqf5 = ddqf{5};  ddqf6 = ddqf{6}; ddqf7 = ddqf{7}; ddqf8 = ddqf{8};
tf1 = tf{1};  tf2 = tf{2}; tf3 = tf{3}; tf4 = tf{4};
tf5 = tf{5};  tf6 = tf{6}; tf7 = tf{7}; tf8 = tf{8};


% % 画出滤波后的加速度曲线
% figure,hold on;
% subplot(421),hold on;plot(dqf1),plot(dq_rm{1}),legend('vel_f', 'vel_r')
% subplot(422),hold on;plot(dqf2),plot(dq_rm{2})
% subplot(423),hold on;plot(dqf3),plot(dq_rm{3})
% subplot(424),hold on;plot(dqf4),plot(dq_rm{4})
% subplot(425),hold on;plot(dqf5),plot(dq_rm{5})
% subplot(426),hold on;plot(dqf6),plot(dq_rm{6})
% subplot(427),hold on;plot(dqf7),plot(dq_rm{7})
% subplot(428),hold on;plot(dqf8),plot(dq_rm{8})
% title('vel')
% 
% 
% % 画出滤波后的速度曲线
% figure,hold on;
% subplot(421),hold on;plot(ddqf1),plot(ddqf1_),legend('acc_f', 'acc_r')
% subplot(422),hold on;plot(ddqf2),plot(ddqf2_)
% subplot(423),hold on;plot(ddqf3),plot(ddqf3_)
% subplot(424),hold on;plot(ddqf4),plot(ddqf4_)
% subplot(425),hold on;plot(ddqf5),plot(ddqf5_)
% subplot(426),hold on;plot(ddqf6),plot(ddqf6_)
% subplot(427),hold on;plot(ddqf7),plot(ddqf7_)
% subplot(428),hold on;plot(ddqf8),plot(ddqf8_)
% title('acc')
% 
% 
% % 画出滤波后的力矩曲线
% figure,hold on;
% subplot(421),hold on;plot(tf1),plot(t_rm{1}),legend('tor_f', 'tor_r')
% subplot(422),hold on;plot(tf2),plot(t_rm{2})
% subplot(423),hold on;plot(tf3),plot(t_rm{3})
% subplot(424),hold on;plot(tf4),plot(t_rm{4})
% subplot(425),hold on;plot(tf5),plot(t_rm{5})
% subplot(426),hold on;plot(tf6),plot(t_rm{6})
% subplot(427),hold on;plot(tf7),plot(t_rm{7})
% subplot(428),hold on;plot(tf8),plot(t_rm{8})
% title('tor')


for n = 1:numOfData
    i = step*(n-1)+si;
    q = [qf1(i), qf2(i),qf3(i),qf4(i),qf5(i),qf6(i),qf7(i),qf8(i)];
    dq = [dqf1(i), dqf2(i),dqf3(i),dqf4(i),dqf5(i),dqf6(i),dqf7(i),dqf8(i)];
    ddq = [ddqf1(i), ddqf2(i),ddqf3(i),ddqf4(i),ddqf5(i),ddqf6(i),ddqf7(i),ddqf8(i)];

    Yr = Pmin_calc(q,dq,ddq);
    Yc = Yc_calc(dq);
    Yt = diag([1 1 1 1 1 1 1 1]);
    Y=[Y;[Yr,Yc,Yt]];

    tau_temp = [tf1(i), tf2(i),tf3(i),tf4(i),tf5(i),tf6(i),tf7(i),tf8(i)];
    tau=[tau, tau_temp];
    tau_real = [tau_real, tau_temp'];
    disp(['proceeding1: ', num2str(n/numOfData*100), ' %'])
end
tau = tau';

P_min_dh=((Y.')*Y)\(Y.')*tau; %最小惯性参数集

tau2 = [];
for n = 1:numOfData
    i = step*(n-1)+si;
    q = [qf1(i), qf2(i),qf3(i),qf4(i),qf5(i),qf6(i),qf7(i),qf8(i)];
    dq = [dqf1(i), dqf2(i),dqf3(i),dqf4(i),dqf5(i),dqf6(i),dqf7(i),dqf8(i)];
    ddq = [ddqf1(i), ddqf2(i),ddqf3(i),ddqf4(i),ddqf5(i),ddqf6(i),ddqf7(i),ddqf8(i)];

    Yr = Pmin_calc(q,dq,ddq);
    Yc = Yc_calc(dq);
    Yt = diag([1 1 1 1 1 1 1 1]);
    tau_tmp = [Yr,Yc,Yt]*P_min_dh;
    tau2=[tau2, tau_tmp];
    disp(['proceeding2: ', num2str(n/numOfData*100), ' %'])
end


figure,title('力矩分开对比图')
subplot(4,4,1),plot(tau2(1,:)),xlabel('index'),ylabel('tor1_iden');
subplot(4,4,2),plot(tau_real(1,:)),xlabel('index'),ylabel('tor1_sim');
subplot(4,4,3),plot(tau2(2,:)),xlabel('index'),ylabel('tor2_iden');
subplot(4,4,4),plot(tau_real(2,:)),xlabel('index'),ylabel('tor2_sim');
subplot(4,4,5),plot(tau2(3,:)),xlabel('index'),ylabel('tor3_iden');
subplot(4,4,6),plot(tau_real(3,:)),xlabel('index'),ylabel('tor3_sim');
subplot(4,4,7),plot(tau2(4,:)),xlabel('index'),ylabel('tor4_iden');
subplot(4,4,8),plot(tau_real(4,:)),xlabel('index'),ylabel('tor4_sim');
subplot(4,4,9),plot(tau2(5,:)),xlabel('index'),ylabel('tor5_iden');
subplot(4,4,10),plot(tau_real(5,:)),xlabel('index'),ylabel('tor5_sim');
subplot(4,4,11),plot(tau2(6,:)),xlabel('index'),ylabel('tor6_iden');
subplot(4,4,12),plot(tau_real(6,:)),xlabel('index'),ylabel('tor6_sim');
subplot(4,4,13),plot(tau2(7,:)),xlabel('index'),ylabel('tor7_iden');
subplot(4,4,14),plot(tau_real(7,:)),xlabel('index'),ylabel('tor7_sim');
subplot(4,4,15),plot(tau2(8,:)),xlabel('index'),ylabel('tor8_iden');
subplot(4,4,16),plot(tau_real(8,:)),xlabel('index'),ylabel('tor8_sim');


figure,title('力矩误差图')
subplot(2,4,1),plot(tau2(1,:)-tau_real(1,:)),xlabel('index'),ylabel('tor_error_1');
subplot(2,4,2),plot(tau2(2,:)-tau_real(2,:)),xlabel('index'),ylabel('tor_error_2');
subplot(2,4,3),plot(tau2(3,:)-tau_real(3,:)),xlabel('index'),ylabel('tor_error_3');
subplot(2,4,4),plot(tau2(4,:)-tau_real(4,:)),xlabel('index'),ylabel('tor_error_4');
subplot(2,4,5),plot(tau2(5,:)-tau_real(5,:)),xlabel('index'),ylabel('tor_error_5');
subplot(2,4,6),plot(tau2(6,:)-tau_real(6,:)),xlabel('index'),ylabel('tor_error_6');
subplot(2,4,7),plot(tau2(7,:)-tau_real(7,:)),xlabel('index'),ylabel('tor_error_7');
subplot(2,4,8),plot(tau2(8,:)-tau_real(8,:)),xlabel('index'),ylabel('tor_error_8');

figure,title('力矩合并对比图')
subplot(2,4,1),hold on;plot(tau2(1,:)),plot(tau_real(1,:)),xlabel('index'),ylabel('tor_1'),legend('td', 'tr')
subplot(2,4,2),hold on;plot(tau2(2,:)),plot(tau_real(2,:)),xlabel('index'),ylabel('tor_2'),legend('td', 'tr')
subplot(2,4,3),hold on;plot(tau2(3,:)),plot(tau_real(3,:)),xlabel('index'),ylabel('tor_3'),legend('td', 'tr')
subplot(2,4,4),hold on;plot(tau2(4,:)),plot(tau_real(4,:)),xlabel('index'),ylabel('tor_4'),legend('td', 'tr')
subplot(2,4,5),hold on;plot(tau2(5,:)),plot(tau_real(5,:)),xlabel('index'),ylabel('tor_5'),legend('td', 'tr')
subplot(2,4,6),hold on;plot(tau2(6,:)),plot(tau_real(6,:)),xlabel('index'),ylabel('tor_6'),legend('td', 'tr')
subplot(2,4,7),hold on;plot(tau2(7,:)),plot(tau_real(7,:)),xlabel('index'),ylabel('tor_7'),legend('td', 'tr')
subplot(2,4,8),hold on;plot(tau2(8,:)),plot(tau_real(8,:)),xlabel('index'),ylabel('tor_8'),legend('td', 'tr')

figure,title('力矩合并对比图')
subplot(2,4,1),hold on;plot(tau2(1,:)),plot(t_rm{1}),xlabel('index'),ylabel('tor_1'),legend('td', 'tr')
subplot(2,4,2),hold on;plot(tau2(2,:)),plot(t_rm{2}),xlabel('index'),ylabel('tor_2'),legend('td', 'tr')
subplot(2,4,3),hold on;plot(tau2(3,:)),plot(t_rm{3}),xlabel('index'),ylabel('tor_3'),legend('td', 'tr')
subplot(2,4,4),hold on;plot(tau2(4,:)),plot(t_rm{4}),xlabel('index'),ylabel('tor_4'),legend('td', 'tr')
subplot(2,4,5),hold on;plot(tau2(5,:)),plot(t_rm{5}),xlabel('index'),ylabel('tor_5'),legend('td', 'tr')
subplot(2,4,6),hold on;plot(tau2(6,:)),plot(t_rm{6}),xlabel('index'),ylabel('tor_6'),legend('td', 'tr')
subplot(2,4,7),hold on;plot(tau2(7,:)),plot(t_rm{7}),xlabel('index'),ylabel('tor_7'),legend('td', 'tr')
subplot(2,4,8),hold on;plot(tau2(8,:)),plot(t_rm{8}),xlabel('index'),ylabel('tor_8'),legend('td', 'tr')


figure
hold on;plot(tau2(1,:)),plot(t_rm{1}),xlabel('index'),ylabel('tor_1'),legend('td', 'tr')
hold on;plot(tau2(2,:)),plot(t_rm{2}),xlabel('index'),ylabel('tor_2'),legend('td', 'tr')
hold on;plot(tau2(3,:)),plot(t_rm{3}),xlabel('index'),ylabel('tor_3'),legend('td', 'tr')
hold on;plot(tau2(4,:)),plot(t_rm{4}),xlabel('index'),ylabel('tor_4'),legend('td', 'tr')
hold on;plot(tau2(5,:)),plot(t_rm{5}),xlabel('index'),ylabel('tor_5'),legend('td', 'tr')
hold on;plot(tau2(6,:)),plot(t_rm{6}),xlabel('index'),ylabel('tor_6'),legend('td', 'tr')
hold on;plot(tau2(7,:)),plot(t_rm{7}),xlabel('index'),ylabel('tor_7'),legend('td', 'tr')
hold on;plot(tau2(8,:)),plot(t_rm{8}),xlabel('index'),ylabel('tor_8'),legend('td', 'tr')


function [qf, dqf, ddqf, tf, q_rm, dq_rm, t_rm] = data_parse(data,dlta_t_msrd)

% 原始关节角度，角速度，电流数据
q_rm = {data(:,1), data(:,2), data(:,3), data(:,4), data(:,5), data(:,6), data(:,7), data(:,8)};
dq_rm = {data(:,9), data(:,10), data(:,11), data(:,12), data(:,13), data(:,14), data(:,15), data(:,16)};
t_rm = {data(:,17), data(:,18), data(:,19), data(:,20), data(:,21), data(:,22), data(:,23), data(:,24)};

% Wrist coupling (motor -> URDF) for j7/j8
% [q7, q8]^T = [0.5 -0.5; 0.5 0.5] * [m7, m8]^T
wrist_A = [0.5, -0.5; 0.5, 0.5];
m7 = q_rm{7}; m8 = q_rm{8};
q_rm{7} = wrist_A(1,1) * m7 + wrist_A(1,2) * m8;
q_rm{8} = wrist_A(2,1) * m7 + wrist_A(2,2) * m8;

dm7 = dq_rm{7}; dm8 = dq_rm{8};
dq_rm{7} = wrist_A(1,1) * dm7 + wrist_A(1,2) * dm8;
dq_rm{8} = wrist_A(2,1) * dm7 + wrist_A(2,2) * dm8;

% If motor torques are logged, map to URDF joint torques:
% tau_urdf = inv(wrist_A') * tau_motor = [1 -1; 1 1] * [t7; t8]
apply_wrist_coupling_to_torque = true;
if apply_wrist_coupling_to_torque
    t7 = t_rm{7}; t8 = t_rm{8};
    t_rm{7} = t7 - t8;
    t_rm{8} = t7 + t8;
end

qf1 = q_rm{1}; qf2 = q_rm{2}; qf3 = q_rm{3};
qf4 = q_rm{4}; qf5 = q_rm{5}; qf6 = q_rm{6};
qf7 = q_rm{7}; qf8 = q_rm{8};


% ---------------------------------------------------------------------
% Filtering Velocities
% ---------------------------------------------------------------------
% Design filter
vel_filt = designfilt('lowpassiir','FilterOrder',5, ...
        'HalfPowerFrequency',0.005,'DesignMethod','butter');
    
dqf1 = filtfilt(vel_filt,dq_rm{1}); dqf2 = filtfilt(vel_filt,dq_rm{2}); dqf3 = filtfilt(vel_filt,dq_rm{3});
dqf4 = filtfilt(vel_filt,dq_rm{4}); dqf5 = filtfilt(vel_filt,dq_rm{5}); dqf6 = filtfilt(vel_filt,dq_rm{6});
dqf7 = filtfilt(vel_filt,dq_rm{7}); dqf8 = filtfilt(vel_filt,dq_rm{8});

% ------------------------------------------------------------------------
% Estimating accelerations
% ------------------------------------------------------------------------
% Three point central difference(三点中心差分法)
ddqf1_ = zeros(size(dq_rm{1})); ddqf2_ = zeros(size(dq_rm{1})); ddqf3_ = zeros(size(dq_rm{1}));
ddqf4_ = zeros(size(dq_rm{1})); ddqf5_ = zeros(size(dq_rm{1})); ddqf6_ = zeros(size(dq_rm{1}));
ddqf7_ = zeros(size(dq_rm{1})); ddqf8_ = zeros(size(dq_rm{1})); 
for i = 2:length(ddqf1_)-1
   dlta_qd_fltrd_1 =  dqf1(i+1) - dqf1(i-1);  dlta_qd_fltrd_2 =  dqf2(i+1) - dqf2(i-1);
   dlta_qd_fltrd_3 =  dqf3(i+1) - dqf3(i-1);  dlta_qd_fltrd_4 =  dqf4(i+1) - dqf4(i-1);
   dlta_qd_fltrd_5 =  dqf5(i+1) - dqf5(i-1);  dlta_qd_fltrd_6 =  dqf6(i+1) - dqf6(i-1);
   dlta_qd_fltrd_7 =  dqf7(i+1) - dqf7(i-1);  dlta_qd_fltrd_8 =  dqf8(i+1) - dqf8(i-1);
   
   ddqf1_(i) = dlta_qd_fltrd_1/dlta_t_msrd;  ddqf2_(i) = dlta_qd_fltrd_2/dlta_t_msrd;
   ddqf3_(i) = dlta_qd_fltrd_3/dlta_t_msrd;  ddqf4_(i) = dlta_qd_fltrd_4/dlta_t_msrd;
   ddqf5_(i) = dlta_qd_fltrd_5/dlta_t_msrd;  ddqf6_(i) = dlta_qd_fltrd_6/dlta_t_msrd;
   ddqf7_(i) = dlta_qd_fltrd_7/dlta_t_msrd;  ddqf8_(i) = dlta_qd_fltrd_8/dlta_t_msrd;
end

% Zeros phase filtering acceleration obtained by finite difference
accel_filt = designfilt('lowpassiir','FilterOrder',5, ...
        'HalfPowerFrequency',0.005,'DesignMethod','butter');

ddqf1 = filtfilt(accel_filt,ddqf1_); ddqf2 = filtfilt(accel_filt,ddqf2_); ddqf3 = filtfilt(accel_filt,ddqf3_);
ddqf4 = filtfilt(accel_filt,ddqf4_); ddqf5 = filtfilt(accel_filt,ddqf5_); ddqf6 = filtfilt(accel_filt,ddqf6_);
ddqf7 = filtfilt(accel_filt,ddqf7_); ddqf8 = filtfilt(accel_filt,ddqf8_);

% ------------------------------------------------------------------------
% Filtering current
% ------------------------------------------------------------------------
% Zeros phase filtering acceleration obtained by finite difference
curr_filt = designfilt('lowpassiir','FilterOrder',5, ...
        'HalfPowerFrequency',0.0055,'DesignMethod','butter');

tf1 = filtfilt(curr_filt,t_rm{1}); tf2 = filtfilt(curr_filt,t_rm{2}); tf3 = filtfilt(curr_filt,t_rm{3});
tf4 = filtfilt(curr_filt,t_rm{4}); tf5 = filtfilt(curr_filt,t_rm{5}); tf6 = filtfilt(curr_filt,t_rm{6});
tf7 = filtfilt(curr_filt,t_rm{7}); tf8 = filtfilt(curr_filt,t_rm{8});

% tf1 = t_rm{1}; tf2 =t_rm{2}; tf3 = t_rm{3};
% tf4 = t_rm{4}; tf5 = t_rm{5}; tf6 = t_rm{6};
% tf7 = t_rm{7}; tf8 = t_rm{8};

qf = {qf1, qf2, qf3, qf4, qf5, qf6, qf7, qf8};
dqf = {dqf1, dqf2, dqf3, dqf4, dqf5, dqf6, dqf7, dqf8};
ddqf = {ddqf1, ddqf2, ddqf3, ddqf4, ddqf5, ddqf6, ddqf7, ddqf8};
tf = {tf1, tf2, tf3, tf4, tf5, tf6, tf7, tf8};

end

        
