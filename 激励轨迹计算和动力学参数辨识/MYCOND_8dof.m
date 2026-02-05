function [cond_val]=MYCOND(x)
%%计算回归矩阵的条件数
n_joints = 8;           % 关节数
T = 10;                 % 周期时间 (s)
dt = 0.05;              % 采样时间间隔
t = 0:dt:T;             % 时间序列
N = 5;                  % 傅里叶级数阶数
freq_base = 1/T;        % 基频 (Hz)

a = zeros(n_joints, N); % 傅里叶系数a
b = zeros(n_joints, N); % 傅里叶系数b
qs = zeros(n_joints, 1);% 初始位置

cond_val=100;


% MDH参数
alpha1 = 1.57080000;
alpha2 = 1.57079633;
alpha3 = 1.57079633 ;
alpha4 = 1.57080000;
alpha5 = 1.57079265 ;
alpha6 = 1.57080000 ;
alpha7 = 1.57079633;
alpha8 = 1.57079265;

a1 = 0.00000000 ;
a2 = 0.17000012 ;
a3 = -0.15799988;
a4 = 0.00036301 ;
a5 = 0.00000010 ;
a6 = -0.05803391;
a7 = 0.00000000 ;
a8 = 0.00000000 ;

d1 = 0.02675000 ;
d2 = 0.00087438 ;
d3 = -0.00087521;
d4 = 0.35330945 ;
d5 = 0.09229371 ;
d6 = -0.22358993;
d7 = 0.00049999 ;
d8 = 0.00000000 ;


theata10 = -0.00000367 ;  
theata20 = -1.57078898 ;
theata30 = 3.13999633    ;
theata40 = 3.14158898    ;
theata50 = 0.00159265    ;
theata60 = -0.00000000;   
theata70 = -1.57078163;
theata80 = 1.57079633 ;


% theta  d   a  aptha  offect
DH = [0, d1, a1, alpha1, theata10;
    0, d2, a2, alpha2, theata20;
    0, d3, a3, alpha3, theata30;
    0, d4, a4, alpha4, theata40;
    0, d5, a5, alpha5, theata50;
    0, d6, a6, alpha6, theata60;
    0, d7, a7, alpha7, theata70;
    0, d8, a8, alpha8, theata80];



for i=1:n_joints
    for j=1:N
        a(i,j)=x(i,j);
        b(i,j)=x(i,j+N);
    end
    qs(i)=x(i,2*N+1);
end


%这里先计算激励轨迹
for j = 1:length(t)
    for i = 1:n_joints
        sum_q = 0;
        sum_dq = 0;
        sum_ddq = 0;
        for k = 1:N
            w_k = 2*pi * k * freq_base; % 角频率
            sum_q = sum_q + a(i,k)/w_k * sin(w_k*t(j)) - b(i,k)/w_k * cos(w_k*t(j));
            sum_dq = sum_dq + a(i,k) * cos(w_k*t(j)) + b(i,k) * sin(w_k*t(j));
            sum_ddq = sum_ddq -a(i,k)*w_k * sin(w_k*t(j)) + b(i,k)*w_k * cos(w_k*t(j));
        end
        q(i,j) = qs(i) + sum_q;
        dq(i,j) = sum_dq;
        ddq(i,j) = sum_ddq;
    end
end


G=[0,0,9.81]';
YY=[];
% 这里在根据激励轨迹计算回归矩阵  这里选取了100个点
if n_joints==2
    DH =    [0  0    0  0 0;
         0  0  0.5  0 0;
         0   0 0.7   0 0];
     for j=1:length(t)
       Yi=Dynamic_Linear(DH,q(:,j),dq(:,j),ddq(:,j),G,2); 
      MYi=[Yi(:,10),Yi(:,12),Yi(:,13),Yi(:,20)];%  这里在加上摩擦力矩阵项
      YY=[YY;MYi]; 
     end
     cond_val =cond(YY)
elseif n_joints==3
else %n_joints==8

     for j=1:length(t)
	 Yr = Pmin_calc(q(:,j),dq(:,j),ddq(:,j));
         YY=[YY;Yr]; 
     end
     cond_val =cond(YY);

end


end


