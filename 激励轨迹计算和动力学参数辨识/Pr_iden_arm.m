            

addpath("dyn_iden_fun/");

Y = [];
tau = [];
P_min_tar = [];

si = 1;
ei = 20001;
step = 1;
numOfData = (ei-si)/step+1;

tau_real = [];
tau_model = [];
tau_model_ = [];

for n = 1:numOfData
    i = step*(n-1)+si;
    q = [getdatasamples(out.q1,i), getdatasamples(out.q2,i), getdatasamples(out.q3,i),getdatasamples(out.q4,i), getdatasamples(out.q5,i), getdatasamples(out.q6,i), getdatasamples(out.q7,i)];
    dq = [getdatasamples(out.dq1,i), getdatasamples(out.dq2,i), getdatasamples(out.dq3,i),getdatasamples(out.dq4,i), getdatasamples(out.dq5,i), getdatasamples(out.dq6,i), getdatasamples(out.dq7,i)];
    ddq = [getdatasamples(out.ddq1,i), getdatasamples(out.ddq2,i), getdatasamples(out.ddq3,i),getdatasamples(out.ddq4,i), getdatasamples(out.ddq5,i), getdatasamples(out.ddq6,i), getdatasamples(out.ddq7,i)];

    Yr = Pmin_calc(q,dq,ddq);
    Yc = Yc_calc(dq);
    Y=[Y;[Yr,Yc]];

    tau_temp = [getdatasamples(out.t1,i),getdatasamples(out.t2,i),getdatasamples(out.t3,i),getdatasamples(out.t4,i),getdatasamples(out.t5,i),getdatasamples(out.t6,i), getdatasamples(out.t7,i)];
    tau=[tau, tau_temp];
    tau_real = [tau_real, tau_temp'];
    disp(['proceeding1: ', num2str(n/numOfData*100), ' %'])
end
tau = tau';

P_min_dh=((Y.')*Y)\(Y.')*tau; %最小惯性参数集

tau2 = [];
for n = 1:numOfData
    i = step*(n-1)+si;
    q = [getdatasamples(out.q1,i), getdatasamples(out.q2,i), getdatasamples(out.q3,i),getdatasamples(out.q4,i), getdatasamples(out.q5,i), getdatasamples(out.q6,i), getdatasamples(out.q7,i)];
    dq = [getdatasamples(out.dq1,i), getdatasamples(out.dq2,i), getdatasamples(out.dq3,i),getdatasamples(out.dq4,i), getdatasamples(out.dq5,i), getdatasamples(out.dq6,i), getdatasamples(out.dq7,i)];
    ddq = [getdatasamples(out.ddq1,i), getdatasamples(out.ddq2,i), getdatasamples(out.ddq3,i),getdatasamples(out.ddq4,i), getdatasamples(out.ddq5,i), getdatasamples(out.ddq6,i), getdatasamples(out.ddq7,i)];
 
    Yr = Pmin_calc(q,dq,ddq);
    Yc = Yc_calc(dq);
    tau_tmp = [Yr,Yc]*P_min_dh;
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


figure,title('力矩误差图')
subplot(2,4,1),plot(tau2(1,:)-tau_real(1,:)),xlabel('index'),ylabel('tor_error_1');
subplot(2,4,2),plot(tau2(2,:)-tau_real(2,:)),xlabel('index'),ylabel('tor_error_2');
subplot(2,4,3),plot(tau2(3,:)-tau_real(3,:)),xlabel('index'),ylabel('tor_error_3');
subplot(2,4,4),plot(tau2(4,:)-tau_real(4,:)),xlabel('index'),ylabel('tor_error_4');
subplot(2,4,5),plot(tau2(5,:)-tau_real(5,:)),xlabel('index'),ylabel('tor_error_5');
subplot(2,4,6),plot(tau2(6,:)-tau_real(6,:)),xlabel('index'),ylabel('tor_error_6');
subplot(2,4,7),plot(tau2(7,:)-tau_real(7,:)),xlabel('index'),ylabel('tor_error_7');


figure,title('力矩合并对比图')
subplot(2,4,1),hold on;plot(tau2(1,:)),plot(tau_real(1,:)),xlabel('index'),ylabel('tor_1'),legend('td', 'tr')
subplot(2,4,2),hold on;plot(tau2(2,:)),plot(tau_real(2,:)),xlabel('index'),ylabel('tor_2'),legend('td', 'tr')
subplot(2,4,3),hold on;plot(tau2(3,:)),plot(tau_real(3,:)),xlabel('index'),ylabel('tor_3'),legend('td', 'tr')
subplot(2,4,4),hold on;plot(tau2(4,:)),plot(tau_real(4,:)),xlabel('index'),ylabel('tor_4'),legend('td', 'tr')
subplot(2,4,5),hold on;plot(tau2(5,:)),plot(tau_real(5,:)),xlabel('index'),ylabel('tor_5'),legend('td', 'tr')
subplot(2,4,6),hold on;plot(tau2(6,:)),plot(tau_real(6,:)),xlabel('index'),ylabel('tor_6'),legend('td', 'tr')
subplot(2,4,7),hold on;plot(tau2(7,:)),plot(tau_real(7,:)),xlabel('index'),ylabel('tor_7'),legend('td', 'tr')

