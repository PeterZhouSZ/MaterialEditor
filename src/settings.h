#ifndef SETTINGS_H
#define SETTINGS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <set>
#include <string>

#define _MODEL_PREFIX_ "../InvMat/models/"

using vec3 = Eigen::Vector3d;
using mat3 = Eigen::Matrix3d;
using DMatrix = Eigen::MatrixXd;
using IMatrix = Eigen::MatrixXi;
using DSparse = Eigen::SparseMatrix<double>;
using DVector = Eigen::VectorXd;
using IVector = Eigen::VectorXi;
using std::cout;
using std::endl;
using std::string;

inline void tetrahedronize(string modelName)
{
  string command = "../InvMat/libs/tetgen/tetgen -pq5.414a.1 ../InvMat/models/" + modelName + "/" + modelName + ".ply";
  system(command.c_str());
}

#endif // SETTINGS_H
