cmake_minimum_required(VERSION 3.12)
project(robot_dynamics)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找Eigen库
find_package(Eigen3 REQUIRED)

# 查找pybind11库
find_package(pybind11 REQUIRED)

# ============================ 第一步：创建核心动态库 ============================
# 创建共享库（动态库），包含核心的robot_dynamics实现
add_library(robot_dynamics_core SHARED
    robot_dynamics.cpp
)

# 为核心库添加Eigen头文件包含路径
target_include_directories(robot_dynamics_core PUBLIC
    ${EIGEN3_INCLUDE_DIRS}
    # 添加当前目录，确保头文件能被其他文件找到
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 链接Eigen（Eigen通常只有头文件，这一步主要是兼容性）
target_link_libraries(robot_dynamics_core PRIVATE
    ${EIGEN3_LIBRARIES}
)

# ============================ 第二步：创建Python绑定模块 ============================
# 添加共享库
pybind11_add_module(robot_dynamics_module
    robot_dynamics.cpp
    robot_wrapper.cpp
)

# 包含目录
target_include_directories(robot_dynamics_module PRIVATE
    ${EIGEN3_INCLUDE_DIRS}
)

# 链接库
target_link_libraries(robot_dynamics_module PRIVATE
    ${EIGEN3_LIBRARIES}
)