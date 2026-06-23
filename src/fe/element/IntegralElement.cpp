#include "fe/element/IntegralElement.h"

void EFEM::IntegralElement::initial(std::shared_ptr<Element> element)
{
    m_element = element;
    for (int i = 0; i < m_element->nodes.size(); i++) {
        auto node = m_element->nodes[i];
        m_infos.push_back(std::static_pointer_cast<ConstraintNodeInfo>(node->extra_info["constraint_info"]));
    }
    int gp_num = getGaussPointNum();
    for (int index = 0; index < gp_num; index++) {
        Eigen::Matrix3d J0 = getGaussPointJacobi(index);
        m_gauss_J0.push_back(J0);
        Eigen::Matrix3d gp_jacobi_inv0 = J0.inverse();
        m_gauss_J0_inv.push_back(gp_jacobi_inv0);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> gp_gradient = getGaussPointShapeFuncLocalGradient(index).transpose();
        m_gauss_T.push_back(gp_jacobi_inv0.transpose() * gp_gradient);
    }
}

Eigen::Vector3d EFEM::IntegralElement::getPointGlobalPosition(Eigen::Vector3d local_position)
{
    Eigen::VectorXd shape_func = getPointShapeFunc(local_position);
    Eigen::Vector3d global_position = Eigen::Vector3d::Zero();
    for (int i = 0; i < m_infos.size(); i++) {
        global_position = shape_func[i] * m_infos[i]->position;
    }
    return global_position;
}

Eigen::Matrix3d EFEM::IntegralElement::getPointJacobi(Eigen::Vector3d local_position)
{
    Eigen::MatrixXd point_shape_gradient = getPointShapeFuncLocalGradient(local_position);
    Eigen::MatrixXd node_position;
    node_position.resize(3, m_infos.size());
    for (int j = 0; j < m_infos.size(); j++) {
        node_position.block<3, 1>(0, j) = m_infos[j]->position;
    }
    return node_position * point_shape_gradient;
}

Eigen::VectorXd EFEM::IntegralElement::getGaussPointShapeFunc(int index)
{
    Eigen::Vector3d local_position = getGaussPointLocalPosition(index);
    return getPointShapeFunc(local_position);
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::IntegralElement::getGaussPointShapeFuncLocalGradient(int index)
{
    Eigen::Vector3d local_position = getGaussPointLocalPosition(index);
    return getPointShapeFuncLocalGradient(local_position);
}

Eigen::Matrix3d EFEM::IntegralElement::getGaussPointJacobiInv0(int index)
{
    return m_gauss_J0_inv[index];
}

Eigen::Matrix<double, 3, -1> EFEM::IntegralElement::getGaussPointT(int index)
{
    return m_gauss_T[index];
}

// calculate dx/d£
Eigen::Matrix3d EFEM::IntegralElement::getGaussPointJacobi(int index)
{
    Eigen::Vector3d local_position = getGaussPointLocalPosition(index);
    return getPointJacobi(local_position);
}

Eigen::Matrix3d EFEM::IntegralElement::getGaussPointF(int index)
{
    // method 1: from
    // Bonet B . Nonlinear Continuum Mechanics for Finite Element Analysis[J]. 2008.
    // Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
    // Eigen::Matrix3d jacobi0 = getGaussPointJacobi0(index);
    // Eigen::Matrix3d jacobi0_invT = jacobi0.inverse().transpose();
    // Eigen::MatrixXd gp_shape_gradient = getGaussPointShapeFuncLocalGradient(index);
    // for (int i = 0; i < 3; i++) {
    //     Eigen::Vector3d F_vector = Eigen::Vector3d::Zero();
    //     for (int j = 0; j < m_infos.size(); j++) {
    //         Eigen::Vector3d dNX = jacobi0_invT * gp_shape_gradient.block<1, 3>(j, 0).transpose();
    //         F_vector += m_infos[j]->position[i] * dNX;
    //     }
    //     F.block<1, 3>(i, 0) += F_vector.transpose();
    // }
    // return F;

    // method 2: simplification of method 1
    // F= dx/dx_0 = dx/d£ * Inv(dx_0/d£)
    return getGaussPointJacobi(index) * getGaussPointJacobiInv0(index);
}