function Yr = Pmin_calc_6dof(q,dq,ddq)

n_joints = 6;

q1=q(1); q2=q(2); q3=q(3); q4=q(4); q5=q(5); q6=q(6);
dq1=dq(1); dq2=dq(2); dq3=dq(3); dq4=dq(4); dq5=dq(5); dq6=dq(6);
ddq1=ddq(1); ddq2=ddq(2); ddq3=ddq(3); ddq4=ddq(4); ddq5=ddq(5); ddq6=ddq(6);

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

% theta  d   a  alpha
DH = [theata10, d1, a1, alpha1;
    theata20, d2, a2, alpha2;
    theata30, d3, a3, alpha3;
    theata40, d4, a4, alpha4;
    theata50, d5, a5, alpha5;
    theata60, d6, a6, alpha6];

g = 9.80665;

axis1=[0 0 1]'; axis2=[0 0 1]'; axis3=[0 0 1]';
axis4=[0 0 1]'; axis5=[0 0 1]'; axis6=[0 0 1]';

T01 = MDHTrans(DH(1,4), DH(1,3), DH(1,2), q1+DH(1,1));
T12 = MDHTrans(DH(2,4), DH(2,3), DH(2,2), q2+DH(2,1));
T23 = MDHTrans(DH(3,4), DH(3,3), DH(3,2), q3+DH(3,1));
T34 = MDHTrans(DH(4,4), DH(4,3), DH(4,2), q4+DH(4,1));
T45 = MDHTrans(DH(5,4), DH(5,3), DH(5,2), q5+DH(5,1));
T56 = MDHTrans(DH(6,4), DH(6,3), DH(6,2), q6+DH(6,1));

R01=T01(1:3,1:3); R10=R01'; P01=T01(1:3,4);
R12=T12(1:3,1:3); R21=R12'; P12=T12(1:3,4);
R23=T23(1:3,1:3); R32=R23'; P23=T23(1:3,4);
R34=T34(1:3,1:3); R43=R34'; P34=T34(1:3,4);
R45=T45(1:3,1:3); R54=R45'; P45=T45(1:3,4);
R56=T56(1:3,1:3); R65=R56'; P56=T56(1:3,4);

% Forward recursion
w00=[0 0 0]'; dw00=[0 0 0]'; dv00=[0 0 g]'; Pc=[0 0 0]';
[w11,dw11,dv11, ~]=motion_para_clc(R10,dq1,ddq1,w00,dw00,dv00,P01,Pc,axis1);
[w22,dw22,dv22, ~]=motion_para_clc(R21,dq2,ddq2,w11,dw11,dv11,P12,Pc,axis2);
[w33,dw33,dv33, ~]=motion_para_clc(R32,dq3,ddq3,w22,dw22,dv22,P23,Pc,axis3);
[w44,dw44,dv44, ~]=motion_para_clc(R43,dq4,ddq4,w33,dw33,dv33,P34,Pc,axis4);
[w55,dw55,dv55, ~]=motion_para_clc(R54,dq5,ddq5,w44,dw44,dv44,P45,Pc,axis5);
[w66,dw66,dv66, ~]=motion_para_clc(R65,dq6,ddq6,w55,dw55,dv55,P56,Pc,axis6);

% Helper matrices
H1=getHi(w11,dw11,dv11); H2=getHi(w22,dw22,dv22); H3=getHi(w33,dw33,dv33);
H4=getHi(w44,dw44,dv44); H5=getHi(w55,dw55,dv55); H6=getHi(w66,dw66,dv66);

A1=getAi(w11,dw11,dv11); A2=getAi(w22,dw22,dv22); A3=getAi(w33,dw33,dv33);
A4=getAi(w44,dw44,dv44); A5=getAi(w55,dw55,dv55); A6=getAi(w66,dw66,dv66);

% Regressor assembly
n_params = n_joints * 10;
Yf6=[zeros(3,n_params-10) H6]; Yn6=[zeros(3,n_params-10) A6];
[Yf5, Yn5]=get_Yf_Yn(5,R56,H5,A5,Yf6,Yn6,P56,n_params);
[Yf4, Yn4]=get_Yf_Yn(4,R45,H4,A4,Yf5,Yn5,P45,n_params);
[Yf3, Yn3]=get_Yf_Yn(3,R34,H3,A3,Yf4,Yn4,P34,n_params);
[Yf2, Yn2]=get_Yf_Yn(2,R23,H2,A2,Yf3,Yn3,P23,n_params);
[Yf1, Yn1]=get_Yf_Yn(1,R12,H1,A1,Yf2,Yn2,P12,n_params);

Y = zeros(n_joints, n_params);
Y(1,1:n_params)=axis1'*Yn1;
Y(2,1:n_params)=axis2'*Yn2;
Y(3,1:n_params)=axis3'*Yn3;
Y(4,1:n_params)=axis4'*Yn4;
Y(5,1:n_params)=axis5'*Yn5;
Y(6,1:n_params)=axis6'*Yn6;

% R2 R1 R1 R1 R1 R1 (first 6 joints)
Yr=zeros(n_joints,38);
Yr(:,1)=Y(:,6); Yr(:,2)=Y(:,7); Yr(:,3)=Y(:,8);
Yr(:,4)=Y(:,11); Yr(:,5)=Y(:,12); Yr(:,6)=Y(:,13); Yr(:,7)=Y(:,15); Yr(:,8)=Y(:,16); Yr(:,9)=Y(:,17); Yr(:,10)=Y(:,18);
Yr(:,11)=Y(:,21); Yr(:,12)=Y(:,22); Yr(:,13)=Y(:,23); Yr(:,14)=Y(:,25); Yr(:,15)=Y(:,26); Yr(:,16)=Y(:,27); Yr(:,17)=Y(:,28);
Yr(:,18)=Y(:,31); Yr(:,19)=Y(:,32); Yr(:,20)=Y(:,33); Yr(:,21)=Y(:,35); Yr(:,22)=Y(:,36); Yr(:,23)=Y(:,37); Yr(:,24)=Y(:,38);
Yr(:,25)=Y(:,41); Yr(:,26)=Y(:,42); Yr(:,27)=Y(:,43); Yr(:,28)=Y(:,45); Yr(:,29)=Y(:,46); Yr(:,30)=Y(:,47); Yr(:,31)=Y(:,48);
Yr(:,32)=Y(:,51); Yr(:,33)=Y(:,52); Yr(:,34)=Y(:,53); Yr(:,35)=Y(:,55); Yr(:,36)=Y(:,56); Yr(:,37)=Y(:,57); Yr(:,38)=Y(:,58);

end

function [Yf, Yn] = get_Yf_Yn(id, R, H, A, Yf_next, Yn_next, Po, n_params)
    c1=(id-1)*10+1; c2=(id-1)*10+10;
    Yf=zeros(3,n_params); Yf(:,c1:c2)=H;
    Yf=Yf+R*Yf_next;
    Yn=zeros(3,n_params); Yn(:,c1:c2)=A;
    Yn=Yn+R*Yn_next+getS(Po)*R*Yf_next;
end
