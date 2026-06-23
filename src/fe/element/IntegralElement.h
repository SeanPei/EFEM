#ifndef EFEM_CONSTRAINT_ELEMENT_INTEGRALELEMENT_H
#define EFEM_CONSTRAINT_ELEMENT_INTEGRALELEMENT_H

#include "fe/Constraint.h"

namespace EFEM {
enum IntegralElementType {
    TETRAHEDRA_3D4,
    TETRAHEDRA_3D4R,
    TETRAHEDRA_3D10,
    TETRAHEDRA_3D10R,
    HEXAHEDRA_3D8,
    HEXAHEDRA_3D8R,
    HEXAHEDRA_3D20,
    HEXAHEDRA_3D20R
};

template <int num>
class TetGaussIntegrationPoints {
public:
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);
};

template <int num>
class HexGaussIntegrationPoints {
public:
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);
};
#include "fe/element/IntegralElement.inl"

class IntegralElement {
public:
    IntegralElementType type;
    virtual Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position) = 0;
    virtual Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position) = 0;
    virtual int getGaussPointNum() = 0;
    virtual double getGaussPointWeight(int index) = 0;
    virtual Eigen::Vector3d getGaussPointLocalPosition(int index) = 0;

    Eigen::Vector3d getPointGlobalPosition(Eigen::Vector3d local_position);
    Eigen::Matrix3d getPointJacobi(Eigen::Vector3d local_position);
    Eigen::VectorXd getGaussPointShapeFunc(int index);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getGaussPointShapeFuncLocalGradient(int index);
    Eigen::Matrix3d getGaussPointJacobi0(int index);
    Eigen::Matrix3d getGaussPointJacobiInv0(int index);
    Eigen::Matrix<double, 3, -1> getGaussPointT(int index);
    Eigen::Matrix3d getGaussPointJacobi(int index);
    Eigen::Matrix3d getGaussPointF(int index);

protected:
    void initial(std::shared_ptr<Element> element);
    std::shared_ptr<Element> m_element;
    std::vector<std::shared_ptr<ConstraintNodeInfo>> m_infos;
    std::vector<Eigen::Matrix3d> m_gauss_J0;
    std::vector<Eigen::Matrix3d> m_gauss_J0_inv;
    std::vector<Eigen::Matrix<double, 3, -1>> m_gauss_T;
};

class TetElement3D4 : public IntegralElement {
public:
    TetElement3D4(std::shared_ptr<Element> element);
    Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position);
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);

private:
    TetGaussIntegrationPoints<1> m_gauss_points;
};

class TetElement3D4R : public TetElement3D4 {
public:
    TetElement3D4R(std::shared_ptr<Element> element);
};

class TetElement3D10 : public IntegralElement {
public:
    TetElement3D10(std::shared_ptr<Element> element);
    Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position);
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);

private:
    TetGaussIntegrationPoints<4> m_gauss_points;
};

class TetElement3D10R : public IntegralElement {
public:
    TetElement3D10R(std::shared_ptr<Element> element);
    Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position);
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);

private:
    TetGaussIntegrationPoints<1> m_gauss_points;
};

class HexElement3D8 : public IntegralElement {
public:
    HexElement3D8(std::shared_ptr<Element> element);
    Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position);
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);

private:
    HexGaussIntegrationPoints<8> m_gauss_points;
};

class HexElement3D8R : public IntegralElement {
public:
    HexElement3D8R(std::shared_ptr<Element> element);
    Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position);
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);

private:
    HexGaussIntegrationPoints<1> m_gauss_points;
};

class HexElement3D20 : public IntegralElement {
public:
    HexElement3D20(std::shared_ptr<Element> element);
    Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position);
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);

private:
    HexGaussIntegrationPoints<27> m_gauss_points;
};

class HexElement3D20R : public IntegralElement {
public:
    HexElement3D20R(std::shared_ptr<Element> element);
    Eigen::VectorXd getPointShapeFunc(Eigen::Vector3d local_position);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getPointShapeFuncLocalGradient(Eigen::Vector3d local_position);
    int getGaussPointNum();
    double getGaussPointWeight(int index);
    Eigen::Vector3d getGaussPointLocalPosition(int index);

private:
    HexGaussIntegrationPoints<8> m_gauss_points;
};
}

#endif