cmake_minimum_required(VERSION 3.12)
project(robot_dynamics_test)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找Eigen库
find_package(Eigen3 REQUIRED)


# 添加测试可执行程序
add_executable(robot_dynamics_test main.cpp)

# 设置包含目录
target_include_directories(robot_dynamics_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIRS}
)

target_link_directories(robot_dynamics_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib                  # 自己编译的robot_dynamics_core.so目录    # 第三方动态库目录
)


# 链接核心动态库
target_link_libraries(robot_dynamics_test PRIVATE
    robot_dynamics_core
)

# 对于使用M_PI等数学常量的情况，需要链接数学库（Linux下）
if(UNIX AND NOT APPLE)
    target_link_libraries(robot_dynamics_test PRIVATE m)
endif()