#include "fe/element/IntegralElement.h"
#include "spdlog/spdlog.h"

EFEM::HexElement3D8::HexElement3D8(std::shared_ptr<Element> element)
{
    type = IntegralElementType::HEXAHEDRA_3D8;
    IntegralElement::initial(element);
}

int EFEM::HexElement3D8::getGaussPointNum()
{
    return m_gauss_points.getGaussPointNum();
}

Eigen::VectorXd EFEM::HexElement3D8::getPointShapeFunc(Eigen::Vector3d local_position)
{
    Eigen::VectorXd shape_func;
    shape_func.resize(8);
    shape_func[0] = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    shape_func[1] = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    shape_func[2] = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    shape_func[3] = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    shape_func[4] = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    shape_func[5] = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    shape_func[6] = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    shape_func[7] = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    return shape_func;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::HexElement3D8::getPointShapeFuncLocalGradient(Eigen::Vector3d local_position)
{
    Eigen::Matrix<double, 8, 3> local_gradient;
    local_gradient(0, 0) = -0.125 * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    local_gradient(0, 1) = -0.125 * (1.0 - local_position[0]) * (1.0 - local_position[2]);
    local_gradient(0, 2) = -0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]);
    local_gradient(1, 0) = 0.125 * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    local_gradient(1, 1) = -0.125 * (1.0 + local_position[0]) * (1.0 - local_position[2]);
    local_gradient(1, 2) = -0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]);
    local_gradient(2, 0) = 0.125 * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    local_gradient(2, 1) = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[2]);
    local_gradient(2, 2) = -0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]);
    local_gradient(3, 0) = -0.125 * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    local_gradient(3, 1) = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[2]);
    local_gradient(3, 2) = -0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]);
    local_gradient(4, 0) = -0.125 * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    local_gradient(4, 1) = -0.125 * (1.0 - local_position[0]) * (1.0 + local_position[2]);
    local_gradient(4, 2) = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]);
    local_gradient(5, 0) = 0.125 * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    local_gradient(5, 1) = -0.125 * (1.0 + local_position[0]) * (1.0 + local_position[2]);
    local_gradient(5, 2) = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]);
    local_gradient(6, 0) = 0.125 * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    local_gradient(6, 1) = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[2]);
    local_gradient(6, 2) = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]);
    local_gradient(7, 0) = -0.125 * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    local_gradient(7, 1) = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[2]);
    local_gradient(7, 2) = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]);
    return local_gradient;
}

double EFEM::HexElement3D8::getGaussPointWeight(int index)
{
    return m_gauss_points.getGaussPointWeight(index);
}

Eigen::Vector3d EFEM::HexElement3D8::getGaussPointLocalPosition(int index)
{
    return m_gauss_points.getGaussPointLocalPosition(index);
}

EFEM::HexElement3D8R::HexElement3D8R(std::shared_ptr<Element> element)
{
    type = IntegralElementType::HEXAHEDRA_3D8R;
    IntegralElement::initial(element);
}

int EFEM::HexElement3D8R::getGaussPointNum()
{
    return m_gauss_points.getGaussPointNum();
}

Eigen::VectorXd EFEM::HexElement3D8R::getPointShapeFunc(Eigen::Vector3d local_position)
{
    Eigen::VectorXd shape_func;
    shape_func.resize(8);
    shape_func[0] = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    shape_func[1] = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    shape_func[2] = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    shape_func[3] = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    shape_func[4] = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    shape_func[5] = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    shape_func[6] = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    shape_func[7] = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    return shape_func;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::HexElement3D8R::getPointShapeFuncLocalGradient(Eigen::Vector3d local_position)
{
    Eigen::Matrix<double, 8, 3> local_gradient;
    local_gradient(0, 0) = -0.125 * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    local_gradient(0, 1) = -0.125 * (1.0 - local_position[0]) * (1.0 - local_position[2]);
    local_gradient(0, 2) = -0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]);
    local_gradient(1, 0) = 0.125 * (1.0 - local_position[1]) * (1.0 - local_position[2]);
    local_gradient(1, 1) = -0.125 * (1.0 + local_position[0]) * (1.0 - local_position[2]);
    local_gradient(1, 2) = -0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]);
    local_gradient(2, 0) = 0.125 * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    local_gradient(2, 1) = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[2]);
    local_gradient(2, 2) = -0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]);
    local_gradient(3, 0) = -0.125 * (1.0 + local_position[1]) * (1.0 - local_position[2]);
    local_gradient(3, 1) = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[2]);
    local_gradient(3, 2) = -0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]);
    local_gradient(4, 0) = -0.125 * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    local_gradient(4, 1) = -0.125 * (1.0 - local_position[0]) * (1.0 + local_position[2]);
    local_gradient(4, 2) = 0.125 * (1.0 - local_position[0]) * (1.0 - local_position[1]);
    local_gradient(5, 0) = 0.125 * (1.0 - local_position[1]) * (1.0 + local_position[2]);
    local_gradient(5, 1) = -0.125 * (1.0 + local_position[0]) * (1.0 + local_position[2]);
    local_gradient(5, 2) = 0.125 * (1.0 + local_position[0]) * (1.0 - local_position[1]);
    local_gradient(6, 0) = 0.125 * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    local_gradient(6, 1) = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[2]);
    local_gradient(6, 2) = 0.125 * (1.0 + local_position[0]) * (1.0 + local_position[1]);
    local_gradient(7, 0) = -0.125 * (1.0 + local_position[1]) * (1.0 + local_position[2]);
    local_gradient(7, 1) = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[2]);
    local_gradient(7, 2) = 0.125 * (1.0 - local_position[0]) * (1.0 + local_position[1]);
    return local_gradient;
}

double EFEM::HexElement3D8R::getGaussPointWeight(int index)
{
    return m_gauss_points.getGaussPointWeight(index);
}

Eigen::Vector3d EFEM::HexElement3D8R::getGaussPointLocalPosition(int index)
{
    return m_gauss_points.getGaussPointLocalPosition(index);
}

EFEM::HexElement3D20::HexElement3D20(std::shared_ptr<Element> element)
{
    type = IntegralElementType::HEXAHEDRA_3D20;
    IntegralElement::initial(element);
}

int EFEM::HexElement3D20::getGaussPointNum()
{
    return m_gauss_points.getGaussPointNum();
}

Eigen::VectorXd EFEM::HexElement3D20::getPointShapeFunc(Eigen::Vector3d local_position)
{
    Eigen::VectorXd shape_func;
    shape_func.resize(20);
    shape_func[0] = -((1.0 + local_position[0]) * (1.0 - local_position[1]) * (2.0 - local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[1] = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * (2.0 - local_position[0] - local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[2] = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]) * (2.0 - local_position[0] - local_position[1] + local_position[2])) / 8.0;
    shape_func[3] = -((1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]) * (2.0 - local_position[0] + local_position[1] + local_position[2])) / 8.0;
    shape_func[4] = -((1.0 - local_position[0]) * (1.0 - local_position[1]) * (2.0 + local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[5] = -((1.0 - local_position[0]) * (1.0 + local_position[1]) * (2.0 + local_position[0] - local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[6] = -((1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]) * (2.0 + local_position[0] - local_position[1] + local_position[2])) / 8.0;
    shape_func[7] = -((1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]) * (2.0 + local_position[0] + local_position[1] + local_position[2])) / 8.0;
    shape_func[8] = ((1.0 + local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[9] = ((1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[10] = ((1.0 + local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 - local_position[2])) / 4.0;
    shape_func[11] = ((1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[12] = ((1.0 - local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[13] = ((1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[14] = ((1.0 - local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 - local_position[2])) / 4.0;
    shape_func[15] = ((1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[16] = ((1.0 - local_position[0] * local_position[0]) * (1.0 - local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[17] = ((1.0 - local_position[0] * local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[18] = ((1.0 - local_position[0] * local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2])) / 4.0;
    shape_func[19] = ((1.0 - local_position[0] * local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2])) / 4.0;
    return shape_func;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::HexElement3D20::getPointShapeFuncLocalGradient(Eigen::Vector3d local_position)
{
    Eigen::Matrix<double, 20, 3> local_gradient;
    local_gradient(0, 0) = ((-1.0 + local_position[1]) * (1.0 - 2.0 * local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(0, 1) = -((1.0 + local_position[0]) * (1.0 + local_position[2]) * (-1.0 + local_position[0] - 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(0, 2) = -((1.0 + local_position[0]) * (-1.0 + local_position[1]) * (-1.0 + local_position[0] - local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(1, 0) = ((1.0 + local_position[1]) * (1.0 + local_position[2]) * (-1.0 + 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(1, 1) = ((1.0 + local_position[0]) * (1.0 + local_position[2]) * (-1.0 + local_position[0] + 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(1, 2) = ((1.0 + local_position[0]) * (1.0 + local_position[1]) * (-1.0 + local_position[0] + local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(2, 0) = -((1.0 + local_position[1]) * (-1.0 + 2.0 * local_position[0] + local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(2, 1) = -((1.0 + local_position[0]) * (-1.0 + local_position[0] + 2.0 * local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(2, 2) = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * (-1.0 + local_position[0] + local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(3, 0) = -((-1.0 + local_position[1]) * (-1.0 + local_position[2]) * (1.0 - 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(3, 1) = ((1.0 + local_position[0]) * (-1.0 + local_position[0] - 2.0 * local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(3, 2) = ((1.0 + local_position[0]) * (-1.0 + local_position[1]) * (-1.0 + local_position[0] - local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(4, 0) = -((-1.0 + local_position[1]) * (1.0 + 2.0 * local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(4, 1) = -((-1.0 + local_position[0]) * (1.0 + local_position[0] + 2.0 * local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(4, 2) = -((-1.0 + local_position[0]) * (-1.0 + local_position[1]) * (1.0 + local_position[0] + local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(5, 0) = -((1.0 + local_position[1]) * (1.0 + local_position[2]) * (-1.0 - 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(5, 1) = ((-1.0 + local_position[0]) * (1.0 + local_position[0] - 2.0 * local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(5, 2) = ((-1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[0] - local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(6, 0) = ((1.0 + local_position[1]) * (-1.0 - 2.0 * local_position[0] + local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(6, 1) = -((-1.0 + local_position[0]) * (-1.0 + local_position[2]) * (1.0 + local_position[0] - 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(6, 2) = -((-1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[0] - local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(7, 0) = ((-1.0 + local_position[1]) * (-1.0 + local_position[2]) * (1.0 + 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(7, 1) = ((-1.0 + local_position[0]) * (-1.0 + local_position[2]) * (1.0 + local_position[0] + 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(7, 2) = ((-1.0 + local_position[0]) * (-1.0 + local_position[1]) * (1.0 + local_position[0] + local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(8, 0) = -((-1.0 + local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(8, 1) = -((1.0 + local_position[0]) * local_position[1] * (1.0 + local_position[2])) / 2.0;
    local_gradient(8, 2) = -((1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(9, 0) = -((1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(9, 1) = -((1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(9, 2) = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(10, 0) = ((-1.0 + local_position[1] * local_position[1]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(10, 1) = ((1.0 + local_position[0]) * local_position[1] * (-1.0 + local_position[2])) / 2.0;
    local_gradient(10, 2) = ((1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(11, 0) = ((-1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(11, 1) = ((1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(11, 2) = ((1.0 + local_position[0]) * (-1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(12, 0) = ((-1.0 + local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(12, 1) = ((-1.0 + local_position[0]) * local_position[1] * (1.0 + local_position[2])) / 2.0;
    local_gradient(12, 2) = ((-1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(13, 0) = ((1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(13, 1) = ((-1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(13, 2) = ((-1.0 + local_position[0]) * (1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(14, 0) = -((-1.0 + local_position[1] * local_position[1]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(14, 1) = -((-1.0 + local_position[0]) * local_position[1] * (-1.0 + local_position[2])) / 2.0;
    local_gradient(14, 2) = -((-1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(15, 0) = -((-1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(15, 1) = -((-1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(15, 2) = -((-1.0 + local_position[0]) * (-1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(16, 0) = (local_position[0] * (-1.0 + local_position[1]) * (1.0 + local_position[2])) / 2.0;
    local_gradient(16, 1) = ((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(16, 2) = ((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[1])) / 4.0;
    local_gradient(17, 0) = -(local_position[0] * (1.0 + local_position[1]) * (1.0 + local_position[2])) / 2.0;
    local_gradient(17, 1) = -((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(17, 2) = -((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[1])) / 4.0;
    local_gradient(18, 0) = (local_position[0] * (1.0 + local_position[1]) * (-1.0 + local_position[2])) / 2.0;
    local_gradient(18, 1) = ((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(18, 2) = ((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[1])) / 4.0;
    local_gradient(19, 0) = -(local_position[0] * (-1.0 + local_position[1]) * (-1.0 + local_position[2])) / 2.0;
    local_gradient(19, 1) = -((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(19, 2) = -((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[1])) / 4.0;
    return local_gradient;
}

double EFEM::HexElement3D20::getGaussPointWeight(int index)
{
    return m_gauss_points.getGaussPointWeight(index);
}

Eigen::Vector3d EFEM::HexElement3D20::getGaussPointLocalPosition(int index)
{
    return m_gauss_points.getGaussPointLocalPosition(index);
}

EFEM::HexElement3D20R::HexElement3D20R(std::shared_ptr<Element> element)
{
    type = IntegralElementType::HEXAHEDRA_3D20;
    IntegralElement::initial(element);
}

int EFEM::HexElement3D20R::getGaussPointNum()
{
    return m_gauss_points.getGaussPointNum();
}

Eigen::VectorXd EFEM::HexElement3D20R::getPointShapeFunc(Eigen::Vector3d local_position)
{
    Eigen::VectorXd shape_func;
    shape_func.resize(20);
    shape_func[0] = -((1.0 + local_position[0]) * (1.0 - local_position[1]) * (2.0 - local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[1] = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * (2.0 - local_position[0] - local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[2] = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]) * (2.0 - local_position[0] - local_position[1] + local_position[2])) / 8.0;
    shape_func[3] = -((1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]) * (2.0 - local_position[0] + local_position[1] + local_position[2])) / 8.0;
    shape_func[4] = -((1.0 - local_position[0]) * (1.0 - local_position[1]) * (2.0 + local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[5] = -((1.0 - local_position[0]) * (1.0 + local_position[1]) * (2.0 + local_position[0] - local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    shape_func[6] = -((1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2]) * (2.0 + local_position[0] - local_position[1] + local_position[2])) / 8.0;
    shape_func[7] = -((1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2]) * (2.0 + local_position[0] + local_position[1] + local_position[2])) / 8.0;
    shape_func[8] = ((1.0 + local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[9] = ((1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[10] = ((1.0 + local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 - local_position[2])) / 4.0;
    shape_func[11] = ((1.0 + local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[12] = ((1.0 - local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[13] = ((1.0 - local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[14] = ((1.0 - local_position[0]) * (1.0 - local_position[1] * local_position[1]) * (1.0 - local_position[2])) / 4.0;
    shape_func[15] = ((1.0 - local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2] * local_position[2])) / 4.0;
    shape_func[16] = ((1.0 - local_position[0] * local_position[0]) * (1.0 - local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[17] = ((1.0 - local_position[0] * local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[2])) / 4.0;
    shape_func[18] = ((1.0 - local_position[0] * local_position[0]) * (1.0 + local_position[1]) * (1.0 - local_position[2])) / 4.0;
    shape_func[19] = ((1.0 - local_position[0] * local_position[0]) * (1.0 - local_position[1]) * (1.0 - local_position[2])) / 4.0;
    return shape_func;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::HexElement3D20R::getPointShapeFuncLocalGradient(Eigen::Vector3d local_position)
{
    Eigen::Matrix<double, 20, 3> local_gradient;
    local_gradient(0, 0) = ((-1.0 + local_position[1]) * (1.0 - 2.0 * local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(0, 1) = -((1.0 + local_position[0]) * (1.0 + local_position[2]) * (-1.0 + local_position[0] - 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(0, 2) = -((1.0 + local_position[0]) * (-1.0 + local_position[1]) * (-1.0 + local_position[0] - local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(1, 0) = ((1.0 + local_position[1]) * (1.0 + local_position[2]) * (-1.0 + 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(1, 1) = ((1.0 + local_position[0]) * (1.0 + local_position[2]) * (-1.0 + local_position[0] + 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(1, 2) = ((1.0 + local_position[0]) * (1.0 + local_position[1]) * (-1.0 + local_position[0] + local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(2, 0) = -((1.0 + local_position[1]) * (-1.0 + 2.0 * local_position[0] + local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(2, 1) = -((1.0 + local_position[0]) * (-1.0 + local_position[0] + 2.0 * local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(2, 2) = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * (-1.0 + local_position[0] + local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(3, 0) = -((-1.0 + local_position[1]) * (-1.0 + local_position[2]) * (1.0 - 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(3, 1) = ((1.0 + local_position[0]) * (-1.0 + local_position[0] - 2.0 * local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(3, 2) = ((1.0 + local_position[0]) * (-1.0 + local_position[1]) * (-1.0 + local_position[0] - local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(4, 0) = -((-1.0 + local_position[1]) * (1.0 + 2.0 * local_position[0] + local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(4, 1) = -((-1.0 + local_position[0]) * (1.0 + local_position[0] + 2.0 * local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(4, 2) = -((-1.0 + local_position[0]) * (-1.0 + local_position[1]) * (1.0 + local_position[0] + local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(5, 0) = -((1.0 + local_position[1]) * (1.0 + local_position[2]) * (-1.0 - 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(5, 1) = ((-1.0 + local_position[0]) * (1.0 + local_position[0] - 2.0 * local_position[1] - local_position[2]) * (1.0 + local_position[2])) / 8.0;
    local_gradient(5, 2) = ((-1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[0] - local_position[1] - 2.0 * local_position[2])) / 8.0;
    local_gradient(6, 0) = ((1.0 + local_position[1]) * (-1.0 - 2.0 * local_position[0] + local_position[1] - local_position[2]) * (-1.0 + local_position[2])) / 8.0;
    local_gradient(6, 1) = -((-1.0 + local_position[0]) * (-1.0 + local_position[2]) * (1.0 + local_position[0] - 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(6, 2) = -((-1.0 + local_position[0]) * (1.0 + local_position[1]) * (1.0 + local_position[0] - local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(7, 0) = ((-1.0 + local_position[1]) * (-1.0 + local_position[2]) * (1.0 + 2.0 * local_position[0] + local_position[1] + local_position[2])) / 8.0;
    local_gradient(7, 1) = ((-1.0 + local_position[0]) * (-1.0 + local_position[2]) * (1.0 + local_position[0] + 2.0 * local_position[1] + local_position[2])) / 8.0;
    local_gradient(7, 2) = ((-1.0 + local_position[0]) * (-1.0 + local_position[1]) * (1.0 + local_position[0] + local_position[1] + 2.0 * local_position[2])) / 8.0;
    local_gradient(8, 0) = -((-1.0 + local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(8, 1) = -((1.0 + local_position[0]) * local_position[1] * (1.0 + local_position[2])) / 2.0;
    local_gradient(8, 2) = -((1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(9, 0) = -((1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(9, 1) = -((1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(9, 2) = -((1.0 + local_position[0]) * (1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(10, 0) = ((-1.0 + local_position[1] * local_position[1]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(10, 1) = ((1.0 + local_position[0]) * local_position[1] * (-1.0 + local_position[2])) / 2.0;
    local_gradient(10, 2) = ((1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(11, 0) = ((-1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(11, 1) = ((1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(11, 2) = ((1.0 + local_position[0]) * (-1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(12, 0) = ((-1.0 + local_position[1] * local_position[1]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(12, 1) = ((-1.0 + local_position[0]) * local_position[1] * (1.0 + local_position[2])) / 2.0;
    local_gradient(12, 2) = ((-1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(13, 0) = ((1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(13, 1) = ((-1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(13, 2) = ((-1.0 + local_position[0]) * (1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(14, 0) = -((-1.0 + local_position[1] * local_position[1]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(14, 1) = -((-1.0 + local_position[0]) * local_position[1] * (-1.0 + local_position[2])) / 2.0;
    local_gradient(14, 2) = -((-1.0 + local_position[0]) * (-1.0 + local_position[1] * local_position[1])) / 4.0;
    local_gradient(15, 0) = -((-1.0 + local_position[1]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(15, 1) = -((-1.0 + local_position[0]) * (-1.0 + local_position[2] * local_position[2])) / 4.0;
    local_gradient(15, 2) = -((-1.0 + local_position[0]) * (-1.0 + local_position[1]) * local_position[2]) / 2.0;
    local_gradient(16, 0) = (local_position[0] * (-1.0 + local_position[1]) * (1.0 + local_position[2])) / 2.0;
    local_gradient(16, 1) = ((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(16, 2) = ((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[1])) / 4.0;
    local_gradient(17, 0) = -(local_position[0] * (1.0 + local_position[1]) * (1.0 + local_position[2])) / 2.0;
    local_gradient(17, 1) = -((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[2])) / 4.0;
    local_gradient(17, 2) = -((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[1])) / 4.0;
    local_gradient(18, 0) = (local_position[0] * (1.0 + local_position[1]) * (-1.0 + local_position[2])) / 2.0;
    local_gradient(18, 1) = ((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(18, 2) = ((-1.0 + local_position[0] * local_position[0]) * (1.0 + local_position[1])) / 4.0;
    local_gradient(19, 0) = -(local_position[0] * (-1.0 + local_position[1]) * (-1.0 + local_position[2])) / 2.0;
    local_gradient(19, 1) = -((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[2])) / 4.0;
    local_gradient(19, 2) = -((-1.0 + local_position[0] * local_position[0]) * (-1.0 + local_position[1])) / 4.0;
    return local_gradient;
}

double EFEM::HexElement3D20R::getGaussPointWeight(int index)
{
    return m_gauss_points.getGaussPointWeight(index);
}

Eigen::Vector3d EFEM::HexElement3D20R::getGaussPointLocalPosition(int index)
{
    return m_gauss_points.getGaussPointLocalPosition(index);
}