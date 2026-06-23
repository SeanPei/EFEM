#ifndef EFEM_MATH_MATHTYPE_H
#define EFEM_MATH_MATHTYPE_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace EFEM {
typedef double ScalarType;
typedef int IndexType;
typedef Eigen::Matrix<ScalarType, 4, 4, 0, 4, 4> EigenMatrix4;
typedef Eigen::Matrix<ScalarType, 4, 1, 0, 4, 1> EigenVector4;
typedef Eigen::Matrix<ScalarType, 3, 3, 0, 3, 3> EigenMatrix3;
typedef Eigen::Matrix<ScalarType, 3, 1, 0, 3, 1> EigenVector3;
typedef Eigen::Matrix<ScalarType, 2, 2, 0, 2, 2> EigenMatrix2;
typedef Eigen::Matrix<ScalarType, 2, 1, 0, 2, 1> EigenVector2;
typedef Eigen::Matrix<ScalarType, -1, 3, 0, -1, 3> EigenMatrixx3;
typedef Eigen::Matrix<ScalarType, 12, 12, 0, 12, 12> EigenMatrix12;
typedef Eigen::Matrix<ScalarType, 12, 1, 0, 12, 1> EigenVector12;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::SparseMatrix<ScalarType> SparseMatrix;
typedef Eigen::Triplet<ScalarType, IndexType> SparseMatrixTriplet;

class Matrix3333 // 3x3 matrix: each element is a 3x3 matrix
{
public:
    Matrix3333();
    Matrix3333(Matrix3333& other);
    Matrix3333(Eigen::Matrix3d* other);
    ~Matrix3333() {}

    void SetZero(); // [0 0 0; 0 0 0; 0 0 0]; 0 = 3x3 zeros
    void SetIdentity(); //[I 0 0; 0 I 0; 0 0 I]; 0 = 3x3 zeros, I = 3x3 identity

    // operators
    EigenMatrix3& operator()(int row, int col);
    Matrix3333 operator+(const Matrix3333& plus);
    Matrix3333 operator-(const Matrix3333& minus);
    Matrix3333 operator*(const EigenMatrix3& multi);
    friend Matrix3333 operator*(const EigenMatrix3& multi1, const Matrix3333& multi2);
    Matrix3333 operator*(ScalarType multi);
    friend Matrix3333 operator*(ScalarType multi1, const Matrix3333& multi2);
    Matrix3333 transpose();
    EigenMatrix3 Contract(const EigenMatrix3& multi); // this operator is commutative
    Matrix3333 Contract(const Matrix3333& multi);

protected:
    EigenMatrix3 mat[3][3];
};

Matrix3333 operator*(const EigenMatrix3& multi1, const Matrix3333& multi2);
Matrix3333 operator*(ScalarType multi1, const Matrix3333& multi2);
}

#endif