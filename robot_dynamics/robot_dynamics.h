#ifndef ROBOT_DYNAMICS_H
#define ROBOT_DYNAMICS_H

#include <vector>
#include <array>
#include <Eigen/Dense>

using namespace Eigen;

class RobotDynamics {
public:
    RobotDynamics();
    
    // 主计算函数
    MatrixXd Pmin_calc(const VectorXd& q, const VectorXd& dq, const VectorXd& ddq);
    Eigen::Matrix<double, 8, 16> Yc_calc(const Eigen::VectorXd& dq);

    
    // 工具函数
    static Matrix3d getS(const Vector3d& vec);
    static MatrixXd getK(const Vector3d& vec);
    static MatrixXd getHi(const Vector3d& w, const Vector3d& dw, const Vector3d& dv);
    static MatrixXd getAi(const Vector3d& w, const Vector3d& dw, const Vector3d& dv);
    static Matrix4d MDHTrans(double alpha, double a, double d, double theta);
    static void motion_para_clc(const Matrix3d& R_inv, double dq, double ddq,
                               const Vector3d& w_pre, const Vector3d& dw_pre, const Vector3d& dv_pre,
                               const Vector3d& Po, const Vector3d& Pc, const Vector3d& axis,
                               Vector3d& w, Vector3d& dw, Vector3d& dv, Vector3d& dvc);
    static void get_Yf_Yn(int id, const Matrix3d& R, const MatrixXd& H, const MatrixXd& A,
                         const MatrixXd& Yf_next, const MatrixXd& Yn_next, const Vector3d& Po,
                         MatrixXd& Yf, MatrixXd& Yn);
    
private:
    // MDH参数
    std::array<double, 8> alpha = {1.57080000, 1.57080000, 1.57079265, 1.57080000,
                                    1.57080000, 1.57080000, 1.57080000, 1.57080000};
    std::array<double, 8> a = {0.00013190, -0.17000000, 0.15005143, 0.00020163,
                               0.00000000, -0.05816200, -0.00000001, 0.00000000};
    std::array<double, 8> d = {0.03320012, 0.00066988, 0.00029644, 0.35330010,
                               -0.09228969, 0.21680999, -0.00150000, 0.00000000};
    std::array<double, 8> theta0 = {3.14159265, 1.56920367, 0.00000001, 0.00000000,
                                    3.14159265, -1.57080367, -1.57080000, 3.14159265};
    
    // 重力加速度
    const double g = 9.80665;
    
    // 旋转轴
    std::array<Vector3d, 8> axes;
    
    // 初始化函数
    void initAxes();
};

#endif // ROBOT_DYNAMICS_H
