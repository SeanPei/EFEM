#include "math/Math.h"
#include "Eigen/SVD"

void EFEM::Math::singularValueDecomp(Eigen::Matrix3d& U, Eigen::Vector3d& SIGMA, Eigen::Matrix3d& V, const Eigen::Matrix3d& A, bool signed_svd)
{
    // Eigen Jacobi SVD
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;
    svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    U = svd.matrixU();
    V = svd.matrixV();
    SIGMA = svd.singularValues();

    if (signed_svd) {
        double detU = U.determinant();
        double detV = V.determinant();
        if (detU < 0) {
            U.block<3, 1>(0, 2) *= -1;
            SIGMA[2] *= -1;
        }
        if (detV < 0) {
            V.block<3, 1>(0, 2) *= -1;
            SIGMA[2] *= -1;
        }
    }
}

void EFEM::Math::kroneckerProduct(Matrix3333& dst, const EigenMatrix3& src1, const EigenMatrix3& src2)
{
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            dst(i, j) = src1(i, j) * src2;
        }
    }
}

Eigen::Matrix4d EFEM::Math::poseMatrix(const std::vector<std::pair<std::string, std::vector<double>>>& pose)
{
    Eigen::Matrix4d pose_matrix = Eigen::Matrix4d::Identity();
    for (auto pose_pair : pose) {
        if (pose_pair.first == "r") {
            Eigen::Matrix4d r_x, r_y, r_z;
            if (pose_pair.second.size() == 3) {
                double t_x = S_PI * pose_pair.second[0] / 180.0;
                double t_y = S_PI * pose_pair.second[1] / 180.0;
                double t_z = S_PI * pose_pair.second[2] / 180.0;
                r_x << 1, 0, 0, 0,
                    0, cos(t_x), -sin(t_x), 0,
                    0, sin(t_x), cos(t_x), 0,
                    0, 0, 0, 1;
                r_y << cos(t_y), 0, sin(t_y), 0,
                    0, 1, 0, 0,
                    -sin(t_y), 0, cos(t_y), 0,
                    0, 0, 0, 1;
                r_z << cos(t_z), -sin(t_z), 0, 0,
                    sin(t_z), cos(t_z), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
                pose_matrix = r_z * r_y * r_x * pose_matrix;
            } else if (pose_pair.second.size() == 4) {
                double u = pose_pair.second[0];
                double v = pose_pair.second[1];
                double w = pose_pair.second[2];
                double t = S_PI * pose_pair.second[3] / 180.8;
                double l = std::sqrt(u * u + v * v + w * w);
                u = u / l;
                v = v / l;
                w = w / l;
                Eigen::Matrix4d r;
                r << u * u + (1 - u * u) * cos(t), u * v * (1 - cos(t)) - w * sin(t), u * w * (1 - cos(t)) + v * sin(t), 0,
                    u * v * (1 - cos(t)) + w * sin(t), v * v + (1 - v * v) * cos(t), v * w * (1 - cos(t)) - u * sin(t), 0,
                    u * w * (1 - cos(t)) - v * sin(t), v * w * (1 - cos(t)) + u * sin(t), w * w + (1 - w * w) * cos(t), 0,
                    0, 0, 0, 1;
                pose_matrix = r * pose_matrix;
            }
        } else if (pose_pair.first == "s") {
            double s_x = pose_pair.second[0];
            double s_y = pose_pair.second[1];
            double s_z = pose_pair.second[2];
            Eigen::Matrix4d s;
            s << s_x, 0, 0, 0,
                0, s_y, 0, 0,
                0, 0, s_z, 0,
                0, 0, 0, 1;
            pose_matrix = s * pose_matrix;
        } else if (pose_pair.first == "t") {
            double t_x = pose_pair.second[0];
            double t_y = pose_pair.second[1];
            double t_z = pose_pair.second[2];
            Eigen::Matrix4d t;
            t << 1, 0, 0, t_x,
                0, 1, 0, t_y,
                0, 0, 1, t_z,
                0, 0, 0, 1;
            pose_matrix = t * pose_matrix;
        }
    }
    return pose_matrix;
}

Eigen::Vector3d EFEM::Math::poseTransform(const std::vector<std::pair<std::string, std::vector<double>>>& pose, const Eigen::Vector3d& position)
{
    Eigen::Matrix4d pose_matrix = poseMatrix(pose);
    Eigen::Vector4d node_pos_homogeneous = Eigen::Vector4d::Ones();
    node_pos_homogeneous.segment(0, 3) = position;
    node_pos_homogeneous = pose_matrix * node_pos_homogeneous;
    return node_pos_homogeneous.segment(0, 3);
}

Eigen::Vector3d EFEM::Math::poseTransform(const Eigen::Matrix4d& pose_matrix, const Eigen::Vector3d& position)
{
    Eigen::Vector4d node_pos_homogeneous = Eigen::Vector4d::Ones();
    node_pos_homogeneous.segment(0, 3) = position;
    node_pos_homogeneous = pose_matrix * node_pos_homogeneous;
    return node_pos_homogeneous.segment(0, 3);
}
