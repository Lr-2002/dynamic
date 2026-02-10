
addpath("dyn_iden_fun/");

% Load measured data
% Update this path as needed
data_path = "traj_1_modify/forward_modify/traj_12-19-18-09-47/uts_traj_left.txt";
data = importdata(data_path);

Y = [];
tau = [];

si = 1;
ei = size(data,1);
step = 1;
numOfData = (ei-si)/step+1;
dlta_t_msrd = 0.002;

tau_real = [];

[qf, dqf, ddqf, tf, q_rm, dq_rm, t_rm] = data_parse(data,dlta_t_msrd);
qf1 = qf{1};  qf2 = qf{2}; qf3 = qf{3}; qf4 = qf{4};
qf5 = qf{5};  qf6 = qf{6};
dqf1 = dqf{1};  dqf2 = dqf{2}; dqf3 = dqf{3}; dqf4 = dqf{4};
dqf5 = dqf{5};  dqf6 = dqf{6};
ddqf1 = ddqf{1};  ddqf2 = ddqf{2}; ddqf3 = ddqf{3}; ddqf4 = ddqf{4};
ddqf5 = ddqf{5};  ddqf6 = ddqf{6};
tf1 = tf{1};  tf2 = tf{2}; tf3 = tf{3}; tf4 = tf{4};
tf5 = tf{5};  tf6 = tf{6};

for n = 1:numOfData
    i = step*(n-1)+si;
    q = [qf1(i), qf2(i), qf3(i), qf4(i), qf5(i), qf6(i)];
    dq = [dqf1(i), dqf2(i), dqf3(i), dqf4(i), dqf5(i), dqf6(i)];
    ddq = [ddqf1(i), ddqf2(i), ddqf3(i), ddqf4(i), ddqf5(i), ddqf6(i)];

    Yr = Pmin_calc_6dof(q,dq,ddq);
    Yc = Yc_calc_6dof(dq);
    Yt = diag(ones(1,6));
    Y = [Y; [Yr, Yc, Yt]];

    tau_temp = [tf1(i), tf2(i), tf3(i), tf4(i), tf5(i), tf6(i)];
    tau = [tau, tau_temp];
    tau_real = [tau_real, tau_temp'];
    disp(['proceeding1: ', num2str(n/numOfData*100), ' %'])
end

tau = tau';

P_min_dh = ((Y.')*Y)\(Y.')*tau;
save('P_min_dh_6dof.mat', 'P_min_dh');

tau2 = [];
for n = 1:numOfData
    i = step*(n-1)+si;
    q = [qf1(i), qf2(i), qf3(i), qf4(i), qf5(i), qf6(i)];
    dq = [dqf1(i), dqf2(i), dqf3(i), dqf4(i), dqf5(i), dqf6(i)];
    ddq = [ddqf1(i), ddqf2(i), ddqf3(i), ddqf4(i), ddqf5(i), ddqf6(i)];

    Yr = Pmin_calc_6dof(q,dq,ddq);
    Yc = Yc_calc_6dof(dq);
    Yt = diag(ones(1,6));
    tau_tmp = [Yr, Yc, Yt] * P_min_dh;
    tau2 = [tau2, tau_tmp];
    disp(['proceeding2: ', num2str(n/numOfData*100), ' %'])
end

figure('Name','Torque compare')
for j = 1:6
    subplot(3,4,2*j-1), plot(tau2(j,:)), xlabel('index'), ylabel(['tor', num2str(j), '_iden']);
    subplot(3,4,2*j), plot(tau_real(j,:)), xlabel('index'), ylabel(['tor', num2str(j), '_real']);
end

figure('Name','Torque error')
for j = 1:6
    subplot(2,3,j), plot(tau2(j,:)-tau_real(j,:)), xlabel('index'), ylabel(['tor_error_', num2str(j)]);
end

figure('Name','Torque compare combined')
for j = 1:6
    subplot(2,3,j), hold on; plot(tau2(j,:)), plot(tau_real(j,:)), xlabel('index'), ylabel(['tor_', num2str(j)]), legend('td', 'tr');
end

figure('Name','Torque compare filtered')
for j = 1:6
    subplot(2,3,j), hold on; plot(tau2(j,:)), plot(t_rm{j}), xlabel('index'), ylabel(['tor_', num2str(j)]), legend('td', 'tr');
end

function [qf, dqf, ddqf, tf, q_rm, dq_rm, t_rm] = data_parse(data,dlta_t_msrd)

% Raw joint angle, velocity, current data
q_rm = {data(:,1), data(:,2), data(:,3), data(:,4), data(:,5), data(:,6), data(:,7), data(:,8)};
dq_rm = {data(:,9), data(:,10), data(:,11), data(:,12), data(:,13), data(:,14), data(:,15), data(:,16)};
t_rm = {data(:,17), data(:,18), data(:,19), data(:,20), data(:,21), data(:,22), data(:,23), data(:,24)};

% Wrist coupling (motor -> URDF) for j7/j8
wrist_A = [0.5, -0.5; 0.5, 0.5];
m7 = q_rm{7}; m8 = q_rm{8};
q_rm{7} = wrist_A(1,1) * m7 + wrist_A(1,2) * m8;
q_rm{8} = wrist_A(2,1) * m7 + wrist_A(2,2) * m8;

dm7 = dq_rm{7}; dm8 = dq_rm{8};
dq_rm{7} = wrist_A(1,1) * dm7 + wrist_A(1,2) * dm8;
dq_rm{8} = wrist_A(2,1) * dm7 + wrist_A(2,2) * dm8;

% Map motor torques to URDF joint torques
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
% Filtering velocities
% ---------------------------------------------------------------------
vel_filt = designfilt('lowpassiir','FilterOrder',5, ...
        'HalfPowerFrequency',0.005,'DesignMethod','butter');



dqf1 = filtfilt(vel_filt,dq_rm{1}); dqf2 = filtfilt(vel_filt,dq_rm{2}); dqf3 = filtfilt(vel_filt,dq_rm{3});
dqf4 = filtfilt(vel_filt,dq_rm{4}); dqf5 = filtfilt(vel_filt,dq_rm{5}); dqf6 = filtfilt(vel_filt,dq_rm{6});
dqf7 = filtfilt(vel_filt,dq_rm{7}); dqf8 = filtfilt(vel_filt,dq_rm{8});

% ------------------------------------------------------------------------
% Estimating accelerations
% ------------------------------------------------------------------------
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

accel_filt = designfilt('lowpassiir','FilterOrder',5, ...
        'HalfPowerFrequency',0.005,'DesignMethod','butter');

ddqf1 = filtfilt(accel_filt,ddqf1_); ddqf2 = filtfilt(accel_filt,ddqf2_); ddqf3 = filtfilt(accel_filt,ddqf3_);
ddqf4 = filtfilt(accel_filt,ddqf4_); ddqf5 = filtfilt(accel_filt,ddqf5_); ddqf6 = filtfilt(accel_filt,ddqf6_);
ddqf7 = filtfilt(accel_filt,ddqf7_); ddqf8 = filtfilt(accel_filt,ddqf8_);

% ------------------------------------------------------------------------
% Filtering current
% ------------------------------------------------------------------------
curr_filt = designfilt('lowpassiir','FilterOrder',5, ...
        'HalfPowerFrequency',0.0055,'DesignMethod','butter');

tf1 = filtfilt(curr_filt,t_rm{1}); tf2 = filtfilt(curr_filt,t_rm{2}); tf3 = filtfilt(curr_filt,t_rm{3});
tf4 = filtfilt(curr_filt,t_rm{4}); tf5 = filtfilt(curr_filt,t_rm{5}); tf6 = filtfilt(curr_filt,t_rm{6});
tf7 = filtfilt(curr_filt,t_rm{7}); tf8 = filtfilt(curr_filt,t_rm{8});

qf = {qf1, qf2, qf3, qf4, qf5, qf6, qf7, qf8};
dqf = {dqf1, dqf2, dqf3, dqf4, dqf5, dqf6, dqf7, dqf8};
ddqf = {ddqf1, ddqf2, ddqf3, ddqf4, ddqf5, ddqf6, ddqf7, ddqf8};
tf = {tf1, tf2, tf3, tf4, tf5, tf6, tf7, tf8};

end
