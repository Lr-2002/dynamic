#include "robot_dynamics.h"
#include <iostream>

RobotDynamics::RobotDynamics() {
    initAxes();
}

void RobotDynamics::initAxes() {
    for (int i = 0; i < 8; i++) {
        axes[i] = Vector3d(0, 0, 1);
    }
}

MatrixXd RobotDynamics::Pmin_calc(const VectorXd& q, const VectorXd& dq, const VectorXd& ddq) {
    // 检查输入维度
    if (q.size() != 8 || dq.size() != 8 || ddq.size() != 8) {
        std::cerr << "错误: 输入向量必须是8维!" << std::endl;
        return MatrixXd::Zero(8, 52);
    }
    
    // 提取关节角度
    double q1 = q(0), q2 = q(1), q3 = q(2), q4 = q(3);
    double q5 = q(4), q6 = q(5), q7 = q(6), q8 = q(7);
    
    // 提取关节速度
    double dq1 = dq(0), dq2 = dq(1), dq3 = dq(2), dq4 = dq(3);
    double dq5 = dq(4), dq6 = dq(5), dq7 = dq(6), dq8 = dq(7);
    
    // 提取关节加速度
    double ddq1 = ddq(0), ddq2 = ddq(1), ddq3 = ddq(2), ddq4 = ddq(3);
    double ddq5 = ddq(4), ddq6 = ddq(5), ddq7 = ddq(6), ddq8 = ddq(7);
    
    // // 打印调试信息
    // std::cout << "q: " << q.transpose() << std::endl;
    // std::cout << "dq: " << dq.transpose() << std::endl;
    // std::cout << "ddq: " << ddq.transpose() << std::endl;
    
    // 计算变换矩阵
    Matrix4d T01 = MDHTrans(alpha[0], a[0], d[0], q1 + theta0[0]);
    Matrix4d T12 = MDHTrans(alpha[1], a[1], d[1], q2 + theta0[1]);
    Matrix4d T23 = MDHTrans(alpha[2], a[2], d[2], q3 + theta0[2]);
    Matrix4d T34 = MDHTrans(alpha[3], a[3], d[3], q4 + theta0[3]);
    Matrix4d T45 = MDHTrans(alpha[4], a[4], d[4], q5 + theta0[4]);
    Matrix4d T56 = MDHTrans(alpha[5], a[5], d[5], q6 + theta0[5]);
    Matrix4d T67 = MDHTrans(alpha[6], a[6], d[6], q7 + theta0[6]);
    Matrix4d T78 = MDHTrans(alpha[7], a[7], d[7], q8 + theta0[7]);
    Matrix4d T89 = MDHTrans(0, 0, 0, 0);
    
    // 提取旋转矩阵和位置向量
    Matrix3d R01 = T01.block<3,3>(0,0);
    Matrix3d R10 = R01.transpose();
    Vector3d P01 = T01.block<3,1>(0,3);
    
    Matrix3d R12 = T12.block<3,3>(0,0);
    Matrix3d R21 = R12.transpose();
    Vector3d P12 = T12.block<3,1>(0,3);
    
    Matrix3d R23 = T23.block<3,3>(0,0);
    Matrix3d R32 = R23.transpose();
    Vector3d P23 = T23.block<3,1>(0,3);
    
    Matrix3d R34 = T34.block<3,3>(0,0);
    Matrix3d R43 = R34.transpose();
    Vector3d P34 = T34.block<3,1>(0,3);
    
    Matrix3d R45 = T45.block<3,3>(0,0);
    Matrix3d R54 = R45.transpose();
    Vector3d P45 = T45.block<3,1>(0,3);
    
    Matrix3d R56 = T56.block<3,3>(0,0);
    Matrix3d R65 = R56.transpose();
    Vector3d P56 = T56.block<3,1>(0,3);
    
    Matrix3d R67 = T67.block<3,3>(0,0);
    Matrix3d R76 = R67.transpose();
    Vector3d P67 = T67.block<3,1>(0,3);
    
    Matrix3d R78 = T78.block<3,3>(0,0);
    Matrix3d R87 = R78.transpose();
    Vector3d P78 = T78.block<3,1>(0,3);
    
    // // 打印调试信息：检查变换矩阵
    // std::cout << "\n调试信息 - 变换矩阵:" << std::endl;
    // std::cout << "R56:\n" << R56 << std::endl;
    // std::cout << "P56: " << P56.transpose() << std::endl;
    
    // 运动参数计算
    Vector3d w00 = Vector3d::Zero();
    Vector3d dw00 = Vector3d::Zero();
    Vector3d dv00(0, 0, g);  // 注意：Matlab代码中dv00是[0 0 g]'
    Vector3d Pc = Vector3d::Zero();
    Vector3d dvc_temp;
    
    Vector3d w11, dw11, dv11;
    motion_para_clc(R10, dq1, ddq1, w00, dw00, dv00, P01, Pc, axes[0], 
                    w11, dw11, dv11, dvc_temp);
    
    Vector3d w22, dw22, dv22;
    motion_para_clc(R21, dq2, ddq2, w11, dw11, dv11, P12, Pc, axes[1],
                    w22, dw22, dv22, dvc_temp);
    
    Vector3d w33, dw33, dv33;
    motion_para_clc(R32, dq3, ddq3, w22, dw22, dv22, P23, Pc, axes[2],
                    w33, dw33, dv33, dvc_temp);
    
    Vector3d w44, dw44, dv44;
    motion_para_clc(R43, dq4, ddq4, w33, dw33, dv33, P34, Pc, axes[3],
                    w44, dw44, dv44, dvc_temp);
    
    Vector3d w55, dw55, dv55;
    motion_para_clc(R54, dq5, ddq5, w44, dw44, dv44, P45, Pc, axes[4],
                    w55, dw55, dv55, dvc_temp);
    
    Vector3d w66, dw66, dv66;
    motion_para_clc(R65, dq6, ddq6, w55, dw55, dv55, P56, Pc, axes[5],
                    w66, dw66, dv66, dvc_temp);
    
    Vector3d w77, dw77, dv77;
    motion_para_clc(R76, dq7, ddq7, w66, dw66, dv66, P67, Pc, axes[6],
                    w77, dw77, dv77, dvc_temp);
    
    Vector3d w88, dw88, dv88;
    motion_para_clc(R87, dq8, ddq8, w77, dw77, dv77, P78, Pc, axes[7],
                    w88, dw88, dv88, dvc_temp);

    // // 打印调试信息：检查运动参数
    // std::cout << "\n调试信息 - 运动参数:" << std::endl;
    // std::cout << "w55: " << w55.transpose() << std::endl;
    // std::cout << "dv55: " << dv55.transpose() << std::endl;
    // std::cout << "dw55: " << dw55.transpose() << std::endl;


    // 计算H和A矩阵
    MatrixXd H1 = getHi(w11, dw11, dv11);
    MatrixXd H2 = getHi(w22, dw22, dv22);
    MatrixXd H3 = getHi(w33, dw33, dv33);
    MatrixXd H4 = getHi(w44, dw44, dv44);
    MatrixXd H5 = getHi(w55, dw55, dv55);
    MatrixXd H6 = getHi(w66, dw66, dv66);
    MatrixXd H7 = getHi(w77, dw77, dv77);
    MatrixXd H8 = getHi(w88, dw88, dv88);
    
    MatrixXd A1 = getAi(w11, dw11, dv11);
    MatrixXd A2 = getAi(w22, dw22, dv22);
    MatrixXd A3 = getAi(w33, dw33, dv33);
    MatrixXd A4 = getAi(w44, dw44, dv44);
    MatrixXd A5 = getAi(w55, dw55, dv55);
    MatrixXd A6 = getAi(w66, dw66, dv66);
    MatrixXd A7 = getAi(w77, dw77, dv77);
    MatrixXd A8 = getAi(w88, dw88, dv88);
    
    // // 打印调试信息：检查H和A矩阵
    // std::cout << "\n调试信息 - H5矩阵形状: " << H1.rows() << "x" << H1.cols() << std::endl;
    // std::cout << "H5:\n" << H5.block<3,10>(0,0) << std::endl;
    // std::cout << "A5:\n" << A5.block<3,10>(0,0) << std::endl;
    
    // 初始化Yf和Yn矩阵 - 修正：应该是80列不是79列
    MatrixXd Yf8 = MatrixXd::Zero(3, 80);
    MatrixXd Yn8 = MatrixXd::Zero(3, 80);
    
    // 注意：Matlab中的列索引是1-based，所以70对应C++中的索引69开始
    // H8是3x10矩阵，应该放到第70-79列（索引69-78）
    Yf8.block<3,10>(0, 70) = H8;
    Yn8.block<3,10>(0, 70) = A8;
    
    // 递归计算Yf和Yn
    MatrixXd Yf7, Yn7;
    get_Yf_Yn(7, R78, H7, A7, Yf8, Yn8, P78, Yf7, Yn7);
    
    MatrixXd Yf6, Yn6;
    get_Yf_Yn(6, R67, H6, A6, Yf7, Yn7, P67, Yf6, Yn6);
    
    MatrixXd Yf5, Yn5;
    get_Yf_Yn(5, R56, H5, A5, Yf6, Yn6, P56, Yf5, Yn5);
    
    MatrixXd Yf4, Yn4;
    get_Yf_Yn(4, R45, H4, A4, Yf5, Yn5, P45, Yf4, Yn4);
    
    MatrixXd Yf3, Yn3;
    get_Yf_Yn(3, R34, H3, A3, Yf4, Yn4, P34, Yf3, Yn3);
    
    MatrixXd Yf2, Yn2;
    get_Yf_Yn(2, R23, H2, A2, Yf3, Yn3, P23, Yf2, Yn2);
    
    MatrixXd Yf1, Yn1;
    get_Yf_Yn(1, R12, H1, A1, Yf2, Yn2, P12, Yf1, Yn1);
    
    // // 打印调试信息：检查Yn1
    // std::cout << "\n调试信息 - Yn5形状: " << Yn5.rows() << "x" << Yn5.cols() << std::endl;
    // std::cout << "Yn8:\n" << Yn8.block<3,10>(0,70) << std::endl;
    // std::cout << "Yf8:\n" << Yf8.block<3,10>(0,70) << std::endl;

    
    // 计算Y矩阵 - 修正：使用正确的转置
    MatrixXd Y = MatrixXd::Zero(8, 80);
    
    // axes是列向量，需要转置后与3x80矩阵相乘，得到1x80行向量
    Y.row(0) = axes[0].transpose() * Yn1;
    Y.row(1) = axes[1].transpose() * Yn2;
    Y.row(2) = axes[2].transpose() * Yn3;
    Y.row(3) = axes[3].transpose() * Yn4;
    Y.row(4) = axes[4].transpose() * Yn5;
    Y.row(5) = axes[5].transpose() * Yn6;
    Y.row(6) = axes[6].transpose() * Yn7;
    Y.row(7) = axes[7].transpose() * Yn8;
    
    // // 打印调试信息：检查Y矩阵
    // std::cout << "\n调试信息 - Y矩阵形状: " << Y.rows() << "x" << Y.cols() << std::endl;
    // std::cout << "Y矩阵第二行前80个元素: " << Y.row(1).head(80) << std::endl;
    
    // 提取Yr矩阵 - 修正索引问题
    // Matlab使用1-based索引，C++使用0-based索引
    MatrixXd Yr = MatrixXd::Zero(8, 52);
    
    // Matlab代码中的索引需要减1转换为C++索引
    // 第1列: Y(:,6) -> Y.col(5)
    // 第2列: Y(:,7) -> Y.col(6)
    // 第3列: Y(:,8) -> Y.col(7)
    Yr.col(0) = Y.col(5);
    Yr.col(1) = Y.col(6);
    Yr.col(2) = Y.col(7);
    
    // 第4-10列
    Yr.col(3) = Y.col(10);   // Y(:,11)
    Yr.col(4) = Y.col(11);   // Y(:,12)
    Yr.col(5) = Y.col(12);   // Y(:,13)
    Yr.col(6) = Y.col(14);   // Y(:,15)
    Yr.col(7) = Y.col(15);   // Y(:,16)
    Yr.col(8) = Y.col(16);   // Y(:,17)
    Yr.col(9) = Y.col(17);   // Y(:,18)
    
    // 第11-17列
    Yr.col(10) = Y.col(20);  // Y(:,21)
    Yr.col(11) = Y.col(21);  // Y(:,22)
    Yr.col(12) = Y.col(22);  // Y(:,23)
    Yr.col(13) = Y.col(24);  // Y(:,25)
    Yr.col(14) = Y.col(25);  // Y(:,26)
    Yr.col(15) = Y.col(26);  // Y(:,27)
    Yr.col(16) = Y.col(27);  // Y(:,28)
    
    // 第18-24列
    Yr.col(17) = Y.col(30);  // Y(:,31)
    Yr.col(18) = Y.col(31);  // Y(:,32)
    Yr.col(19) = Y.col(32);  // Y(:,33)
    Yr.col(20) = Y.col(34);  // Y(:,35)
    Yr.col(21) = Y.col(35);  // Y(:,36)
    Yr.col(22) = Y.col(36);  // Y(:,37)
    Yr.col(23) = Y.col(37);  // Y(:,38)
    
    // 第25-31列
    Yr.col(24) = Y.col(40);  // Y(:,41)
    Yr.col(25) = Y.col(41);  // Y(:,42)
    Yr.col(26) = Y.col(42);  // Y(:,43)
    Yr.col(27) = Y.col(44);  // Y(:,45)
    Yr.col(28) = Y.col(45);  // Y(:,46)
    Yr.col(29) = Y.col(46);  // Y(:,47)
    Yr.col(30) = Y.col(47);  // Y(:,48)
    
    // 第32-38列
    Yr.col(31) = Y.col(50);  // Y(:,51)
    Yr.col(32) = Y.col(51);  // Y(:,52)
    Yr.col(33) = Y.col(52);  // Y(:,53)
    Yr.col(34) = Y.col(54);  // Y(:,55)
    Yr.col(35) = Y.col(55);  // Y(:,56)
    Yr.col(36) = Y.col(56);  // Y(:,57)
    Yr.col(37) = Y.col(57);  // Y(:,58)
    
    // 第39-45列
    Yr.col(38) = Y.col(60);  // Y(:,61)
    Yr.col(39) = Y.col(61);  // Y(:,62)
    Yr.col(40) = Y.col(62);  // Y(:,63)
    Yr.col(41) = Y.col(64);  // Y(:,65)
    Yr.col(42) = Y.col(65);  // Y(:,66)
    Yr.col(43) = Y.col(66);  // Y(:,67)
    Yr.col(44) = Y.col(67);  // Y(:,68)
    
    // 第46-52列
    Yr.col(45) = Y.col(70);  // Y(:,71)
    Yr.col(46) = Y.col(71);  // Y(:,72)
    Yr.col(47) = Y.col(72);  // Y(:,73)
    Yr.col(48) = Y.col(74);  // Y(:,75)
    Yr.col(49) = Y.col(75);  // Y(:,76)
    Yr.col(50) = Y.col(76);  // Y(:,77)
    Yr.col(51) = Y.col(77);  // Y(:,78)
    
    // std::cout << "\nYr矩阵计算完成!" << std::endl;
    // std::cout << "Yr形状: " << Yr.rows() << "x" << Yr.cols() << std::endl;
    // std::cout << "Yr第一行前5个元素: " << Yr.row(0).head(5) << std::endl;
    
    return Yr;
}

// 工具函数实现
Matrix3d RobotDynamics::getS(const Vector3d& vec) {
    Matrix3d S;
    S << 0, -vec(2), vec(1),
         vec(2), 0, -vec(0),
         -vec(1), vec(0), 0;
    return S;
}

MatrixXd RobotDynamics::getK(const Vector3d& vec) {
    MatrixXd K(3, 6);
    K << vec(0), vec(1), vec(2), 0, 0, 0,
         0, vec(0), 0, vec(1), vec(2), 0,
         0, 0, vec(0), 0, vec(1), vec(2);
    return K;
}

MatrixXd RobotDynamics::getHi(const Vector3d& w, const Vector3d& dw, const Vector3d& dv) {
    MatrixXd Hi = MatrixXd::Zero(3, 10);
    
    // 注意：Matlab中是Hi(1:3,7:9) = getS(dwi)+getS(wi)*getS(wi);
    // 对应C++中是第7-9列（索引6-8）
    Hi.block<3,3>(0, 6) = getS(dw) + getS(w) * getS(w);
    
    // Matlab中是Hi(1:3,10) = dvi;
    // 对应C++中是第10列（索引9）
    Hi.col(9) = dv;
    
    return Hi;
}

MatrixXd RobotDynamics::getAi(const Vector3d& w, const Vector3d& dw, const Vector3d& dv) {
    MatrixXd Ai = MatrixXd::Zero(3, 10);
    
    // Matlab中是Ai(1:3,1:6) = getK(dwi)+getS(wi)*getK(wi);
    // 对应C++中是第1-6列（索引0-5）
    Ai.block<3,6>(0, 0) = getK(dw) + getS(w) * getK(w);
    
    // Matlab中是Ai(1:3,7:9) = -getS(dvi);
    // 对应C++中是第7-9列（索引6-8）
    Ai.block<3,3>(0, 6) = -getS(dv);
    
    // Matlab中是Ai(1:3,10) = zeros(3,1);
    // 对应C++中是第10列（索引9）
    Ai.col(9) = Vector3d::Zero();
    
    return Ai;
}

Matrix4d RobotDynamics::MDHTrans(double alpha, double a, double d, double theta) {
    Matrix4d T;
    
    double ct = cos(theta);
    double st = sin(theta);
    double ca = cos(alpha);
    double sa = sin(alpha);
    
    T << ct, -st, 0, a,
         st*ca, ct*ca, -sa, -d*sa,
         st*sa, ct*sa, ca, d*ca,
         0, 0, 0, 1;
    
    return T;
}

void RobotDynamics::motion_para_clc(const Matrix3d& R_inv, double dq, double ddq,
                                   const Vector3d& w_pre, const Vector3d& dw_pre, const Vector3d& dv_pre,
                                   const Vector3d& Po, const Vector3d& Pc, const Vector3d& axis,
                                   Vector3d& w, Vector3d& dw, Vector3d& dv, Vector3d& dvc) {
    // w = R_inv*w_pre+dq*axis;
    w = R_inv * w_pre + dq * axis;
    
    // dw = R_inv*dw_pre+cross(R_inv*w_pre,dq*axis)+ddq*axis;
    dw = R_inv * dw_pre + (R_inv * w_pre).cross(dq * axis) + ddq * axis;
    
    // dv = R_inv*(cross(dw_pre,Po)+cross(w_pre,cross(w_pre,Po))+dv_pre);
    Vector3d temp = dw_pre.cross(Po) + w_pre.cross(w_pre.cross(Po)) + dv_pre;
    dv = R_inv * temp;
    
    // dvc = cross(dw,Pc)+cross(w,cross(w,Pc))+dv;
    dvc = dw.cross(Pc) + w.cross(w.cross(Pc)) + dv;
}

void RobotDynamics::get_Yf_Yn(int id, const Matrix3d& R, const MatrixXd& H, const MatrixXd& A,
                             const MatrixXd& Yf_next, const MatrixXd& Yn_next, const Vector3d& Po,
                             MatrixXd& Yf, MatrixXd& Yn) {
    // 注意：Matlab中id是从1开始的，对应关节编号
    // c1=(id-1)*10+1; c2=(id-1)*10+10;
    // 在C++中，索引从0开始，所以：
    int c1 = (id - 1) * 10;      // 对应Matlab的c1-1
    int c2 = (id - 1) * 10 + 9;             // 对应Matlab的c2-1
    
    Yf = MatrixXd::Zero(3, 80);
    Yf.block<3,10>(0, c1) = H;
    Yf = Yf + R * Yf_next;
    
    Yn = MatrixXd::Zero(3, 80);
    Yn.block<3,10>(0, c1) = A;
    Yn = Yn + R * Yn_next + getS(Po) * R * Yf_next;

    // std::cout << "get_Yf_Yn: " << id << ", " << Yf_next << std::endl;
}

Eigen::Matrix<double, 8, 16> RobotDynamics::Yc_calc(const Eigen::VectorXd& dq) {

    // 定义阈值常量
    const double DQ_MIN = 0.001;

    // 初始化8x16矩阵，所有元素默认置0
    Eigen::Matrix<double, 8, 16> Yc = Eigen::Matrix<double, 8, 16>::Zero();

    // 遍历8个关节维度，逐个计算Kc和Kv并填充矩阵
    for (int i = 0; i < 8; ++i) {
        double dq_i = dq(i);  // 获取第i个关节的角速度
        double Kc = 0.0, Kv = 0.0;

        // 判断角速度绝对值是否大于阈值
        if (std::fabs(dq_i) >= DQ_MIN) {
            // 等价于MATLAB的sign函数：正数=1，负数=-1，0=0
            Kc = (dq_i > 0) ? 1.0 : (dq_i < 0) ? -1.0 : 0.0;
            Kv = dq_i;  // 阈值内则Kv等于角速度本身
        }

        // 填充矩阵：第i行的第2*i列存Kc，第2*i+1列存Kv（与MATLAB逻辑完全一致）
        Yc(i, 2 * i)     = Kc;
        Yc(i, 2 * i + 1) = Kv;
    }

    return Yc;
}
