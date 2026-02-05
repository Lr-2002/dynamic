function Yr = Pmin_calc_func(q,dq,ddq)


q1=q(1); q2=q(2); q3=q(3);q4=q(4); q5=q(5); q6=q(6); q7=q(7); q8=q(8);
dq1=dq(1); dq2=dq(2); dq3=dq(3);dq4=dq(4); dq5=dq(5); dq6=dq(6); dq7=dq(7); dq8=dq(8);
ddq1=ddq(1); ddq2=ddq(2); ddq3=ddq(3);ddq4=ddq(4); ddq5=ddq(5); ddq6=ddq(6); ddq7=ddq(7); ddq8=ddq(8);

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


%  theta      d       a        alpha 
DH = [theata10, d1, a1, alpha1;
    theata20, d2, a2, alpha2;
    theata30, d3, a3, alpha3;
    theata40, d4, a4, alpha4;
    theata50, d5, a5, alpha5;
    theata60, d6, a6, alpha6;
    theata70, d7, a7, alpha7;
    theata80, d8, a8, alpha8];

g = 9.80665;

% 坐标系旋转轴
axis1=[0 0 1]';axis2=[0 0 1]';axis3=[0 0 1]';
axis4=[0 0 1]';axis5=[0 0 1]';axis6=[0 0 1]';axis7=[0 0 1]';axis8=[0 0 1]';

% 坐标系变换矩阵
T01 = MDHTrans(DH(1,4), DH(1,3), DH(1,2), q1+DH(1,1));
T12 = MDHTrans(DH(2,4), DH(2,3), DH(2,2), q2+DH(2,1));
T23 = MDHTrans(DH(3,4), DH(3,3), DH(3,2), q3+DH(3,1));
T34 = MDHTrans(DH(4,4), DH(4,3), DH(4,2), q4+DH(4,1));
T45 = MDHTrans(DH(5,4), DH(5,3), DH(5,2), q5+DH(5,1));
T56 = MDHTrans(DH(6,4), DH(6,3), DH(6,2), q6+DH(6,1));
T67 = MDHTrans(DH(7,4), DH(7,3), DH(7,2), q7+DH(7,1));
T78 = MDHTrans(DH(8,4), DH(8,3), DH(8,2), q8+DH(8,1));
T89 = MDHTrans(0, 0, 0, 0);

R01=T01(1:3,1:3);R10=R01';P01=T01(1:3,4);
R12=T12(1:3,1:3);R21=R12';P12=T12(1:3,4);
R23=T23(1:3,1:3);R32=R23';P23=T23(1:3,4);
R34=T34(1:3,1:3);R43=R34';P34=T34(1:3,4);
R45=T45(1:3,1:3);R54=R45';P45=T45(1:3,4);
R56=T56(1:3,1:3);R65=R56';P56=T56(1:3,4);
R67=T67(1:3,1:3);R76=R67';P67=T67(1:3,4);
R78=T78(1:3,1:3);R87=R78';P78=T78(1:3,4);
R89=T89(1:3,1:3);R98=R89';P89=T89(1:3,4);

%% 外推(dv00里的g的方向和重力的方向相反)
w00=[0 0 0]';dw00=[0 0 0]';dv00=[0 0 g]'; Pc =[0 0 0]';   
[w11,dw11,dv11, ~]=motion_para_clc(R10,dq1,ddq1,w00,dw00,dv00,P01,Pc,axis1);
[w22,dw22,dv22, ~]=motion_para_clc(R21,dq2,ddq2,w11,dw11,dv11,P12,Pc,axis2);
[w33,dw33,dv33, ~]=motion_para_clc(R32,dq3,ddq3,w22,dw22,dv22,P23,Pc,axis3);
[w44,dw44,dv44, ~]=motion_para_clc(R43,dq4,ddq4,w33,dw33,dv33,P34,Pc,axis4);
[w55,dw55,dv55, ~]=motion_para_clc(R54,dq5,ddq5,w44,dw44,dv44,P45,Pc,axis5);
[w66,dw66,dv66, ~]=motion_para_clc(R65,dq6,ddq6,w55,dw55,dv55,P56,Pc,axis6);
[w77,dw77,dv77, ~]=motion_para_clc(R76,dq7,ddq7,w66,dw66,dv66,P67,Pc,axis7);
[w88,dw88,dv88, ~]=motion_para_clc(R87,dq8,ddq8,w77,dw77,dv77,P78,Pc,axis8);

% 辅助矩阵H
H1=getHi(w11,dw11,dv11);
H2=getHi(w22,dw22,dv22);
H3=getHi(w33,dw33,dv33);
H4=getHi(w44,dw44,dv44);
H5=getHi(w55,dw55,dv55);
H6=getHi(w66,dw66,dv66);
H7=getHi(w77,dw77,dv77);
H8=getHi(w88,dw88,dv88);

% 辅助矩阵A
A1=getAi(w11,dw11,dv11);
A2=getAi(w22,dw22,dv22);
A3=getAi(w33,dw33,dv33);
A4=getAi(w44,dw44,dv44);
A5=getAi(w55,dw55,dv55);
A6=getAi(w66,dw66,dv66);
A7=getAi(w77,dw77,dv77);
A8=getAi(w88,dw88,dv88);

% 辅助矩阵Yf、Yn
Yf8=[zeros(3,70) H8]; Yn8=[zeros(3,70) A8]; 
[Yf7, Yn7]=get_Yf_Yn(7,R78,H7,A7,Yf8,Yn8,P78);
[Yf6, Yn6]=get_Yf_Yn(6,R67,H6,A6,Yf7,Yn7,P67);
[Yf5, Yn5]=get_Yf_Yn(5,R56,H5,A5,Yf6,Yn6,P56);
[Yf4, Yn4]=get_Yf_Yn(4,R45,H4,A4,Yf5,Yn5,P45);
[Yf3, Yn3]=get_Yf_Yn(3,R34,H3,A3,Yf4,Yn4,P34);
[Yf2, Yn2]=get_Yf_Yn(2,R23,H2,A2,Yf3,Yn3,P23);
[Yf1, Yn1]=get_Yf_Yn(1,R12,H1,A1,Yf2,Yn2,P12);

% 线性矩阵Y
Y(1,1:80)=axis1'*Yn1;
Y(2,1:80)=axis2'*Yn2;
Y(3,1:80)=axis3'*Yn3;
Y(4,1:80)=axis4'*Yn4;
Y(5,1:80)=axis5'*Yn5;
Y(6,1:80)=axis6'*Yn6;
Y(7,1:80)=axis7'*Yn7;
Y(8,1:80)=axis8'*Yn8;

% R1 ——> 0 ~ (i-1)关节中存在 “和关节 i 轴向不平行” 的转动关节
% Yr=zeros(7,49);
% Yr(:,1)=Y(:,1); Yr(:,2)=Y(:,2); Yr(:,3)=Y(:,3); Yr(:,4)=Y(:,5); Yr(:,5)=Y(:,6); Yr(:,6)=Y(:,7); Yr(:,7)=Y(:,8);
% Yr(:,8)=Y(:,11);Yr(:,9)=Y(:,12);Yr(:,10)=Y(:,13);Yr(:,11)=Y(:,15);Yr(:,12)=Y(:,16);Yr(:,13)=Y(:,17);Yr(:,14)=Y(:,18);
% Yr(:,15)=Y(:,21);Yr(:,16)=Y(:,22);Yr(:,17)=Y(:,23);Yr(:,18)=Y(:,25);Yr(:,19)=Y(:,26);Yr(:,20)=Y(:,27);Yr(:,21)=Y(:,28);
% Yr(:,22)=Y(:,31);Yr(:,23)=Y(:,32);Yr(:,24)=Y(:,33);Yr(:,25)=Y(:,35);Yr(:,26)=Y(:,36);Yr(:,27)=Y(:,37);Yr(:,28)=Y(:,38);
% Yr(:,29)=Y(:,41);Yr(:,30)=Y(:,42);Yr(:,31)=Y(:,43);Yr(:,32)=Y(:,45);Yr(:,33)=Y(:,46);Yr(:,34)=Y(:,47);Yr(:,35)=Y(:,48);
% Yr(:,36)=Y(:,51);Yr(:,37)=Y(:,52);Yr(:,38)=Y(:,53);Yr(:,39)=Y(:,55);Yr(:,40)=Y(:,56);Yr(:,41)=Y(:,57);Yr(:,42)=Y(:,58);
% Yr(:,43)=Y(:,61);Yr(:,44)=Y(:,62);Yr(:,45)=Y(:,63);Yr(:,46)=Y(:,65);Yr(:,47)=Y(:,66);Yr(:,48)=Y(:,67);Yr(:,49)=Y(:,68);

%R2 ——> 0 ~ (i-1)关节中没有 "和关节 i 轴向不平行" 的转动关节，但存在 “和关节 i 轴向平行但不重合” 的转动关节，或者 “和关节 i 轴向不平行” 的移动关节
% Yr=zeros(7,45);
% Yr(:,1)=Y(:,6); Yr(:,2)=Y(:,7); Yr(:,3)=Y(:,8);
% Yr(:,4)=Y(:,11);Yr(:,5)=Y(:,12);Yr(:,6)=Y(:,13);Yr(:,7)=Y(:,15);Yr(:,8)=Y(:,16);Yr(:,9)=Y(:,17);Yr(:,10)=Y(:,18);
% Yr(:,11)=Y(:,21);Yr(:,12)=Y(:,22);Yr(:,13)=Y(:,23);Yr(:,14)=Y(:,25);Yr(:,15)=Y(:,26);Yr(:,16)=Y(:,27);Yr(:,17)=Y(:,28);
% Yr(:,18)=Y(:,31);Yr(:,19)=Y(:,32);Yr(:,20)=Y(:,33);Yr(:,21)=Y(:,35);Yr(:,22)=Y(:,36);Yr(:,23)=Y(:,37);Yr(:,24)=Y(:,38);
% Yr(:,25)=Y(:,41);Yr(:,26)=Y(:,42);Yr(:,27)=Y(:,43);Yr(:,28)=Y(:,45);Yr(:,29)=Y(:,46);Yr(:,30)=Y(:,47);Yr(:,31)=Y(:,48);
% Yr(:,32)=Y(:,51);Yr(:,33)=Y(:,52);Yr(:,34)=Y(:,53);Yr(:,35)=Y(:,55);Yr(:,36)=Y(:,56);Yr(:,37)=Y(:,57);Yr(:,38)=Y(:,58);
% Yr(:,39)=Y(:,61);Yr(:,40)=Y(:,62);Yr(:,41)=Y(:,63);Yr(:,42)=Y(:,65);Yr(:,43)=Y(:,66);Yr(:,44)=Y(:,67);Yr(:,45)=Y(:,68);

% R3 ——> 0 ~ (i-1)关节中，既没有 “和关节 i 轴向平行但不重合” 的转动关节，也没有 “和关节 i 轴向不平行” 的移动关节
% Yr=zeros(7,43);
% Yr(:,1)=Y(:,6); 
% Yr(:,2)=Y(:,11);Yr(:,3)=Y(:,12);Yr(:,4)=Y(:,13);Yr(:,5)=Y(:,15);Yr(:,6)=Y(:,16);Yr(:,7)=Y(:,17);Yr(:,8)=Y(:,18);
% Yr(:,9)=Y(:,21);Yr(:,10)=Y(:,22);Yr(:,11)=Y(:,23);Yr(:,12)=Y(:,25);Yr(:,13)=Y(:,26);Yr(:,14)=Y(:,27);Yr(:,15)=Y(:,28);
% Yr(:,16)=Y(:,31);Yr(:,17)=Y(:,32);Yr(:,18)=Y(:,33);Yr(:,19)=Y(:,35);Yr(:,20)=Y(:,36);Yr(:,21)=Y(:,37);Yr(:,22)=Y(:,38);
% Yr(:,23)=Y(:,41);Yr(:,24)=Y(:,42);Yr(:,25)=Y(:,43);Yr(:,26)=Y(:,45);Yr(:,27)=Y(:,46);Yr(:,28)=Y(:,47);Yr(:,29)=Y(:,48);
% Yr(:,30)=Y(:,51);Yr(:,31)=Y(:,52);Yr(:,32)=Y(:,53);Yr(:,33)=Y(:,55);Yr(:,34)=Y(:,56);Yr(:,35)=Y(:,57);Yr(:,36)=Y(:,58);
% Yr(:,37)=Y(:,61);Yr(:,38)=Y(:,62);Yr(:,39)=Y(:,63);Yr(:,40)=Y(:,65);Yr(:,41)=Y(:,66);Yr(:,42)=Y(:,67);Yr(:,43)=Y(:,68);

% R3 R1 R1 R1 R1 R1 R1 R1
% Yr=zeros(8,50);
% Yr(:,1)=Y(:,6); 
% Yr(:,2)=Y(:,11);Yr(:,3)=Y(:,12);Yr(:,4)=Y(:,13);Yr(:,5)=Y(:,15);Yr(:,6)=Y(:,16);Yr(:,7)=Y(:,17);Yr(:,8)=Y(:,18);
% Yr(:,9)=Y(:,21);Yr(:,10)=Y(:,22);Yr(:,11)=Y(:,23);Yr(:,12)=Y(:,25);Yr(:,13)=Y(:,26);Yr(:,14)=Y(:,27);Yr(:,15)=Y(:,28);
% Yr(:,16)=Y(:,31);Yr(:,17)=Y(:,32);Yr(:,18)=Y(:,33);Yr(:,19)=Y(:,35);Yr(:,20)=Y(:,36);Yr(:,21)=Y(:,37);Yr(:,22)=Y(:,38);
% Yr(:,23)=Y(:,41);Yr(:,24)=Y(:,42);Yr(:,25)=Y(:,43);Yr(:,26)=Y(:,45);Yr(:,27)=Y(:,46);Yr(:,28)=Y(:,47);Yr(:,29)=Y(:,48);
% Yr(:,30)=Y(:,51);Yr(:,31)=Y(:,52);Yr(:,32)=Y(:,53);Yr(:,33)=Y(:,55);Yr(:,34)=Y(:,56);Yr(:,35)=Y(:,57);Yr(:,36)=Y(:,58);
% Yr(:,37)=Y(:,61);Yr(:,38)=Y(:,62);Yr(:,39)=Y(:,63);Yr(:,40)=Y(:,65);Yr(:,41)=Y(:,66);Yr(:,42)=Y(:,67);Yr(:,43)=Y(:,68);
% Yr(:,44)=Y(:,71);Yr(:,45)=Y(:,72);Yr(:,46)=Y(:,73);Yr(:,47)=Y(:,75);Yr(:,48)=Y(:,76);Yr(:,49)=Y(:,77);Yr(:,50)=Y(:,78);

% R2 R1 R1 R1 R1 R1 R1 R1
Yr=zeros(8,52);
Yr(:,1)=Y(:,6); Yr(:,2)=Y(:,7); Yr(:,3)=Y(:,8); 
Yr(:,4)=Y(:,11);Yr(:,5)=Y(:,12);Yr(:,6)=Y(:,13);Yr(:,7)=Y(:,15);Yr(:,8)=Y(:,16);Yr(:,9)=Y(:,17);Yr(:,10)=Y(:,18);
Yr(:,11)=Y(:,21);Yr(:,12)=Y(:,22);Yr(:,13)=Y(:,23);Yr(:,14)=Y(:,25);Yr(:,15)=Y(:,26);Yr(:,16)=Y(:,27);Yr(:,17)=Y(:,28);
Yr(:,18)=Y(:,31);Yr(:,19)=Y(:,32);Yr(:,20)=Y(:,33);Yr(:,21)=Y(:,35);Yr(:,22)=Y(:,36);Yr(:,23)=Y(:,37);Yr(:,24)=Y(:,38);
Yr(:,25)=Y(:,41);Yr(:,26)=Y(:,42);Yr(:,27)=Y(:,43);Yr(:,28)=Y(:,45);Yr(:,29)=Y(:,46);Yr(:,30)=Y(:,47);Yr(:,31)=Y(:,48);
Yr(:,32)=Y(:,51);Yr(:,33)=Y(:,52);Yr(:,34)=Y(:,53);Yr(:,35)=Y(:,55);Yr(:,36)=Y(:,56);Yr(:,37)=Y(:,57);Yr(:,38)=Y(:,58);
Yr(:,39)=Y(:,61);Yr(:,40)=Y(:,62);Yr(:,41)=Y(:,63);Yr(:,42)=Y(:,65);Yr(:,43)=Y(:,66);Yr(:,44)=Y(:,67);Yr(:,45)=Y(:,68);
Yr(:,46)=Y(:,71);Yr(:,47)=Y(:,72);Yr(:,48)=Y(:,73);Yr(:,49)=Y(:,75);Yr(:,50)=Y(:,76);Yr(:,51)=Y(:,77);Yr(:,52)=Y(:,78);


end


function [Yf, Yn] = get_Yf_Yn(id, R, H, A, Yf_next, Yn_next, Po)
    c1=(id-1)*10+1;c2=(id-1)*10+10;
    Yf=zeros(3,80);Yf(:,c1:c2)=H;
    Yf=Yf+R*Yf_next;
    Yn=zeros(3,80);Yn(:,c1:c2)=A;
    Yn=Yn+R*Yn_next+getS(Po)*R*Yf_next;
end


function Ai=getAi(wi,dwi,dvi)
    Ai(1:3,1:6)=getK(dwi)+getS(wi)*getK(wi);
    Ai(1:3,7:9)=-getS(dvi);
    Ai(1:3,10)=zeros(3,1);
end


function Hi=getHi(wi,dwi,dvi)
    Hi(1:3,7:9)=getS(dwi)+getS(wi)*getS(wi);
    Hi(1:3,10)=dvi;
    Hi(1:3,1:6)=zeros(3,6);
end

function I = getInertiaMatrix(Ic_vec)

    Ixx = Ic_vec(1); Iyy = Ic_vec(2); Izz = Ic_vec(3);
    Iyz = Ic_vec(4); Ixz = Ic_vec(5); Ixy = Ic_vec(6);

    I = [Ixx, Ixy, Ixz;
        Ixy, Iyy, Iyz;
        Ixz, Iyz, Izz];

end


function K=getK(vec3)
    K=[vec3(1)  vec3(2)  vec3(3)    0        0        0;
         0      vec3(1)    0      vec3(2)  vec3(3)    0;
         0        0      vec3(1)    0      vec3(2)  vec3(3)];
end


function S=getS(vec3)
    S=[   0     -vec3(3)  vec3(2);
       vec3(3)     0     -vec3(1);
       -vec3(2) vec3(1)     0    ];
end

function T = MDHTrans(alpha, a, d, theta)
    T = [cos(theta)              -sin(theta)             0              a;
         sin(theta)*cos(alpha)  cos(theta)*cos(alpha)   -sin(alpha)    -d*sin(alpha);
         sin(theta)*sin(alpha)  cos(theta)*sin(alpha)   cos(alpha)     d*cos(alpha);
         0                      0                       0              1];
end


function [w, dw ,dv ,dvc] = motion_para_clc(R_inv,dq,ddq,w_pre,dw_pre,dv_pre,Po,Pc,axis)
    w = R_inv*w_pre+dq*axis;
    dw = R_inv*dw_pre+cross(R_inv*w_pre,dq*axis)+ddq*axis;
    dv = R_inv*(cross(dw_pre,Po)+cross(w_pre,cross(w_pre,Po))+dv_pre);
    dvc = cross(dw,Pc)+cross(w,cross(w,Pc))+dv;
end

function R = RotX(theta)

    ct = cos(theta);
    st = sin(theta);

    R = [
        1   0    0
        0   ct  -st
        0   st   ct
        ];

end

function R = RotY(theta)

    ct = cos(theta);
    st = sin(theta);

    R = [
        ct  0   st
        0   1   0
       -st  0   ct
       ];

end

function R = RotZ(theta)

    ct = cos(theta);
    st = sin(theta);

    R = [
        ct  -st  0
        st   ct  0
        0    0   1
        ];

end
