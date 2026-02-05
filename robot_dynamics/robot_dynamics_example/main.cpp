#include "robot_dynamics.h"
#include <iostream>

int main() {
    // ========== 1. 创建RobotDynamics类实例 ==========
    RobotDynamics robot;
    std::cout << "✅ RobotDynamics 对象创建成功" << std::endl;

    // ========== 2. 准备测试输入数据 ==========
    // 8关节机器人的关节角度、角速度、角加速度示例值
    Eigen::VectorXd q(8);     // 关节角度
    Eigen::VectorXd dq(8);    // 关节角速度
    Eigen::VectorXd ddq(8);   // 关节角加速度
    
    // 初始化测试值（可以根据实际需求修改）
    q << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8;
    dq << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08;
    ddq << 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008;

    // ========== 3. 调用成员函数 ==========
    // 计算Pmin
    Eigen::MatrixXd Pmin = robot.Pmin_calc(q, dq, ddq);
    std::cout << "\n📊 Pmin 计算结果 (维度: " << Pmin.rows() << "x" << Pmin.cols() << "):" << std::endl;
    std::cout << Pmin.block(0, 0, 5, 5) << std::endl;  // 只打印前5x5，避免输出过长

    // 计算Yc
    Eigen::Matrix<double, 8, 16> Yc = robot.Yc_calc(dq);
    std::cout << "\n📊 Yc 计算结果 (维度: " << Yc.rows() << "x" << Yc.cols() << "):" << std::endl;
    std::cout << Yc.block(0, 0, 5, 5) << std::endl;

    // ========== 4. 调用静态工具函数 ==========
    // 测试getS函数（反对称矩阵）
    Eigen::Vector3d vec(1.0, 2.0, 3.0);
    Eigen::Matrix3d S = RobotDynamics::getS(vec);
    std::cout << "\n🔧 getS 静态函数测试 (反对称矩阵):" << std::endl;
    std::cout << S << std::endl;

    // 测试MDHTrans函数（MDH变换矩阵）
    double alpha = M_PI / 2.0;  // 90度
    double a = 0.17;            // 连杆长度
    double d = 0.02675;         // 连杆偏距
    double theta = 0.0;         // 关节角度
    Eigen::Matrix4d T = RobotDynamics::MDHTrans(alpha, a, d, theta);
    std::cout << "\n🔧 MDHTrans 静态函数测试 (变换矩阵):" << std::endl;
    std::cout << T << std::endl;

    // ========== 5. 测试复杂静态函数（motion_para_clc） ==========
    Eigen::Matrix3d R_inv = Eigen::Matrix3d::Identity();  // 单位矩阵
    double dq_test = 0.1;
    double ddq_test = 0.01;
    Eigen::Vector3d w_pre(0, 0, 0);
    Eigen::Vector3d dw_pre(0, 0, 0);
    Eigen::Vector3d dv_pre(0, 0, 0);
    Eigen::Vector3d Po(0.1, 0.2, 0.3);
    Eigen::Vector3d Pc(0.05, 0.05, 0.05);
    Eigen::Vector3d axis(0, 0, 1);  // Z轴
    
    // 输出参数
    Eigen::Vector3d w, dw, dv, dvc;
    
    // 调用运动参数计算函数
    RobotDynamics::motion_para_clc(R_inv, dq_test, ddq_test,
                                   w_pre, dw_pre, dv_pre,
                                   Po, Pc, axis,
                                   w, dw, dv, dvc);
    
    std::cout << "\n🔧 motion_para_clc 静态函数测试:" << std::endl;
    std::cout << "w = " << w.transpose() << std::endl;
    std::cout << "dw = " << dw.transpose() << std::endl;
    std::cout << "dv = " << dv.transpose() << std::endl;
    std::cout << "dvc = " << dvc.transpose() << std::endl;

    std::cout << "\n🎉 所有测试函数调用完成！" << std::endl;

    return 0;
}