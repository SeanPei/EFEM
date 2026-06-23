#include "fe/element/IntegralElement.h"
#include "math/Math.h"

EFEM::TetElement3D4::TetElement3D4(std::shared_ptr<Element> element)
{
    type = IntegralElementType::TETRAHEDRA_3D4;
    IntegralElement::initial(element);
}

int EFEM::TetElement3D4::getGaussPointNum()
{
    return m_gauss_points.getGaussPointNum();
}

Eigen::VectorXd EFEM::TetElement3D4::getPointShapeFunc(Eigen::Vector3d local_position)
{
    Eigen::Vector4d ret;
    ret << local_position[0],
        local_position[1],
        local_position[2],
        1 - local_position[0] - local_position[1] - local_position[2];
    return ret;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::TetElement3D4::getPointShapeFuncLocalGradient(Eigen::Vector3d local_position)
{
    Eigen::Matrix<double, 4, 3> local_gradient;
    local_gradient << 1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        -1, -1, -1;
    return local_gradient;
}

double EFEM::TetElement3D4::getGaussPointWeight(int index)
{
    return m_gauss_points.getGaussPointWeight(index);
}

Eigen::Vector3d EFEM::TetElement3D4::getGaussPointLocalPosition(int index)
{
    return m_gauss_points.getGaussPointLocalPosition(index);
}

EFEM::TetElement3D4R::TetElement3D4R(std::shared_ptr<Element> element)
    : TetElement3D4(element)
{
    type = IntegralElementType::TETRAHEDRA_3D4R;
}

EFEM::TetElement3D10::TetElement3D10(std::shared_ptr<Element> element)
{
    type = IntegralElementType::TETRAHEDRA_3D10;
    IntegralElement::initial(element);
}

int EFEM::TetElement3D10::getGaussPointNum()
{
    return m_gauss_points.getGaussPointNum();
}

Eigen::VectorXd EFEM::TetElement3D10::getPointShapeFunc(Eigen::Vector3d local_position)
{
    Eigen::Vector3d rCoordinates = local_position;
    const double fourthCoord = 1.0 - rCoordinates[0] - rCoordinates[1] - rCoordinates[2];
    Eigen::VectorXd shape_func;
    shape_func.resize(10);
    shape_func[0] = fourthCoord * (2.0 * fourthCoord - 1.0);
    shape_func[1] = rCoordinates[0] * (2.0 * rCoordinates[0] - 1.0);
    shape_func[2] = rCoordinates[1] * (2.0 * rCoordinates[1] - 1.0);
    shape_func[3] = rCoordinates[2] * (2.0 * rCoordinates[2] - 1.0);
    shape_func[4] = 4.0 * fourthCoord * rCoordinates[0];
    shape_func[5] = 4.0 * rCoordinates[0] * rCoordinates[1];
    shape_func[6] = 4.0 * rCoordinates[1] * fourthCoord;
    shape_func[7] = 4.0 * rCoordinates[2] * fourthCoord;
    shape_func[8] = 4.0 * rCoordinates[0] * rCoordinates[2];
    shape_func[9] = 4.0 * rCoordinates[1] * rCoordinates[2];
    return shape_func;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::TetElement3D10::getPointShapeFuncLocalGradient(Eigen::Vector3d local_position)
{
    Eigen::Matrix<double, 10, 3> local_gradient;

    double fourthCoord = 1.0 - (local_position[0] + local_position[1] + local_position[2]);
    double fourthCoord_DX = -1.0;
    double fourthCoord_DY = -1.0;
    double fourthCoord_DZ = -1.0;

    local_gradient(0, 0) = (4 * fourthCoord - 1.0) * fourthCoord_DX;
    local_gradient(0, 1) = (4 * fourthCoord - 1.0) * fourthCoord_DY;
    local_gradient(0, 2) = (4 * fourthCoord - 1.0) * fourthCoord_DZ;
    local_gradient(1, 0) = 4 * local_position[0] - 1.0;
    local_gradient(1, 1) = 0.0;
    local_gradient(1, 2) = 0.0;
    local_gradient(2, 0) = 0.0;
    local_gradient(2, 1) = 4 * local_position[1] - 1.0;
    local_gradient(2, 2) = 0.0;
    local_gradient(3, 0) = 0.0;
    local_gradient(3, 1) = 0.0;
    local_gradient(3, 2) = 4 * local_position[2] - 1.0;
    local_gradient(4, 0) = 4 * fourthCoord_DX * local_position[0] + 4 * fourthCoord;
    local_gradient(4, 1) = 4 * fourthCoord_DY * local_position[0];
    local_gradient(4, 2) = 4 * fourthCoord_DZ * local_position[0];
    local_gradient(5, 0) = 4 * local_position[1];
    local_gradient(5, 1) = 4 * local_position[0];
    local_gradient(5, 2) = 0.0;
    local_gradient(6, 0) = 4 * fourthCoord_DX * local_position[1];
    local_gradient(6, 1) = 4 * fourthCoord_DY * local_position[1] + 4 * fourthCoord;
    local_gradient(6, 2) = 4 * fourthCoord_DZ * local_position[1];
    local_gradient(7, 0) = 4 * fourthCoord_DX * local_position[2];
    local_gradient(7, 1) = 4 * fourthCoord_DY * local_position[2];
    local_gradient(7, 2) = 4 * fourthCoord_DZ * local_position[2] + 4 * fourthCoord;
    local_gradient(8, 0) = 4 * local_position[2];
    local_gradient(8, 1) = 0.0;
    local_gradient(8, 2) = 4 * local_position[0];
    local_gradient(9, 0) = 0.0;
    local_gradient(9, 1) = 4 * local_position[2];
    local_gradient(9, 2) = 4 * local_position[1];

    return local_gradient;
}

double EFEM::TetElement3D10::getGaussPointWeight(int index)
{
    return m_gauss_points.getGaussPointWeight(index);
}

Eigen::Vector3d EFEM::TetElement3D10::getGaussPointLocalPosition(int index)
{
    return m_gauss_points.getGaussPointLocalPosition(index);
}

EFEM::TetElement3D10R::TetElement3D10R(std::shared_ptr<Element> element)
{
    type = IntegralElementType::TETRAHEDRA_3D10R;
    IntegralElement::initial(element);
}

int EFEM::TetElement3D10R::getGaussPointNum()
{
    return m_gauss_points.getGaussPointNum();
}

Eigen::VectorXd EFEM::TetElement3D10R::getPointShapeFunc(Eigen::Vector3d local_position)
{
    Eigen::Vector3d rCoordinates = local_position;
    const double fourthCoord = 1.0 - rCoordinates[0] - rCoordinates[1] - rCoordinates[2];
    Eigen::VectorXd shape_func;
    shape_func.resize(10);
    shape_func[0] = fourthCoord * (2.0 * fourthCoord - 1.0);
    shape_func[1] = rCoordinates[0] * (2.0 * rCoordinates[0] - 1.0);
    shape_func[2] = rCoordinates[1] * (2.0 * rCoordinates[1] - 1.0);
    shape_func[3] = rCoordinates[2] * (2.0 * rCoordinates[2] - 1.0);
    shape_func[4] = 4.0 * fourthCoord * rCoordinates[0];
    shape_func[5] = 4.0 * rCoordinates[0] * rCoordinates[1];
    shape_func[6] = 4.0 * rCoordinates[1] * fourthCoord;
    shape_func[7] = 4.0 * rCoordinates[2] * fourthCoord;
    shape_func[8] = 4.0 * rCoordinates[0] * rCoordinates[2];
    shape_func[9] = 4.0 * rCoordinates[1] * rCoordinates[2];
    return shape_func;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EFEM::TetElement3D10R::getPointShapeFuncLocalGradient(Eigen::Vector3d local_position)
{
    Eigen::Matrix<double, 10, 3> local_gradient;

    double fourthCoord = 1.0 - (local_position[0] + local_position[1] + local_position[2]);
    double fourthCoord_DX = -1.0;
    double fourthCoord_DY = -1.0;
    double fourthCoord_DZ = -1.0;

    local_gradient(0, 0) = (4 * fourthCoord - 1.0) * fourthCoord_DX;
    local_gradient(0, 1) = (4 * fourthCoord - 1.0) * fourthCoord_DY;
    local_gradient(0, 2) = (4 * fourthCoord - 1.0) * fourthCoord_DZ;
    local_gradient(1, 0) = 4 * local_position[0] - 1.0;
    local_gradient(1, 1) = 0.0;
    local_gradient(1, 2) = 0.0;
    local_gradient(2, 0) = 0.0;
    local_gradient(2, 1) = 4 * local_position[1] - 1.0;
    local_gradient(2, 2) = 0.0;
    local_gradient(3, 0) = 0.0;
    local_gradient(3, 1) = 0.0;
    local_gradient(3, 2) = 4 * local_position[2] - 1.0;
    local_gradient(4, 0) = 4 * fourthCoord_DX * local_position[0] + 4 * fourthCoord;
    local_gradient(4, 1) = 4 * fourthCoord_DY * local_position[0];
    local_gradient(4, 2) = 4 * fourthCoord_DZ * local_position[0];
    local_gradient(5, 0) = 4 * local_position[1];
    local_gradient(5, 1) = 4 * local_position[0];
    local_gradient(5, 2) = 0.0;
    local_gradient(6, 0) = 4 * fourthCoord_DX * local_position[1];
    local_gradient(6, 1) = 4 * fourthCoord_DY * local_position[1] + 4 * fourthCoord;
    local_gradient(6, 2) = 4 * fourthCoord_DZ * local_position[1];
    local_gradient(7, 0) = 4 * fourthCoord_DX * local_position[2];
    local_gradient(7, 1) = 4 * fourthCoord_DY * local_position[2];
    local_gradient(7, 2) = 4 * fourthCoord_DZ * local_position[2] + 4 * fourthCoord;
    local_gradient(8, 0) = 4 * local_position[2];
    local_gradient(8, 1) = 0.0;
    local_gradient(8, 2) = 4 * local_position[0];
    local_gradient(9, 0) = 0.0;
    local_gradient(9, 1) = 4 * local_position[2];
    local_gradient(9, 2) = 4 * local_position[1];

    return local_gradient;
}

double EFEM::TetElement3D10R::getGaussPointWeight(int index)
{
    return m_gauss_points.getGaussPointWeight(index);
}

Eigen::Vector3d EFEM::TetElement3D10R::getGaussPointLocalPosition(int index)
{
    return m_gauss_points.getGaussPointLocalPosition(index);
}