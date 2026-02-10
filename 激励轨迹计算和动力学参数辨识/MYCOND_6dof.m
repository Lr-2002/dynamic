function [cond_val]=MYCOND_6dof(x)

n_joints = 6;
T = 10;
dt = 0.05;
t = 0:dt:T;
N = 5;
freq_base = 1/T;

a = zeros(n_joints, N);
b = zeros(n_joints, N);
qs = zeros(n_joints, 1);

cond_val = 100;

% MDH parameters (first 6 joints)
alpha1 = 1.57080000;
alpha2 = 1.57080000;
alpha3 = 1.57079265;
alpha4 = 1.57080000;
alpha5 = 1.57080000;
alpha6 = 1.57080000;

a1 = 0.00013190 ;
a2 = -0.17000000 ;
a3 = 0.15005143;
a4 = 0.00020163 ;
a5 = 0.00000000 ;
a6 = -0.05816200;

d1 = 0.03320012 ;
d2 = 0.00066988 ;
d3 = 0.00029644;
d4 = 0.35330010 ;
d5 = -0.09228969 ;
d6 = 0.21680999;

theata10 = 3.14159265 ;  
theata20 = 1.56920367 ;
theata30 = 0.00000001;
theata40 = 0.00000000;
theata50 = 3.14159265;
theata60 = -1.57080367;

% theta  d   a  alpha  offset
DH = [0, d1, a1, alpha1, theata10;
    0, d2, a2, alpha2, theata20;
    0, d3, a3, alpha3, theata30;
    0, d4, a4, alpha4, theata40;
    0, d5, a5, alpha5, theata50;
    0, d6, a6, alpha6, theata60];

for i=1:n_joints
    for j=1:N
        a(i,j)=x(i,j);
        b(i,j)=x(i,j+N);
    end
    qs(i)=x(i,2*N+1);
end

q = zeros(n_joints, length(t));
dq = zeros(n_joints, length(t));
ddq = zeros(n_joints, length(t));

for j = 1:length(t)
    for i = 1:n_joints
        sum_q = 0;
        sum_dq = 0;
        sum_ddq = 0;
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

YY = [];
for j=1:length(t)
    Yr = Pmin_calc_6dof(q(:,j),dq(:,j),ddq(:,j));
    YY = [YY; Yr];
end
cond_val = cond(YY);

end
