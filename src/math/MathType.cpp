#include "math/MathType.h"

// matrix3333
EFEM::Matrix3333::Matrix3333()
{
}

EFEM::Matrix3333::Matrix3333(Matrix3333& other)
{
    for (unsigned int row = 0; row != 3; ++row) {
        for (unsigned int col = 0; col != 3; ++col) {
            mat[row][col] = other.mat[row][col];
        }
    }
}

EFEM::Matrix3333::Matrix3333(Eigen::Matrix3d* other)
{
    for (unsigned int row = 0; row != 3; ++row) {
        for (unsigned int col = 0; col != 3; ++col) {
            mat[row][col] = other[row * 3 + col];
        }
    }
}

void EFEM::Matrix3333::SetZero()
{
    for (unsigned int row = 0; row != 3; ++row) {
        for (unsigned int col = 0; col != 3; ++col) {
            mat[row][col] = EigenMatrix3::Zero();
        }
    }
}

void EFEM::Matrix3333::SetIdentity()
{
    for (unsigned int row = 0; row != 3; ++row) {
        for (unsigned int col = 0; col != 3; ++col) {
            mat[row][col] = EigenMatrix3::Zero();
            mat[row][col](row, col) = 1.0;
        }
    }
}

EFEM::EigenMatrix3& EFEM::Matrix3333::operator()(int row, int col)
{
    assert(row >= 0 && row < 3 && col >= 0 && col < 3);
    return mat[row][col];
}

EFEM::Matrix3333 EFEM::Matrix3333::operator+(const Matrix3333& plus)
{
    Matrix3333 res;
    for (unsigned int row = 0; row != 3; ++row) {
        for (unsigned int col = 0; col != 3; ++col) {
            res.mat[row][col] = mat[row][col] + plus.mat[row][col];
        }
    }
    return res;
}

EFEM::Matrix3333 EFEM::Matrix3333::operator-(const Matrix3333& minus)
{
    Matrix3333 res;
    for (unsigned int row = 0; row != 3; ++row) {
        for (unsigned int col = 0; col != 3; ++col) {
            res.mat[row][col] = mat[row][col] - minus.mat[row][col];
        }
    }
    return res;
}

EFEM::Matrix3333 EFEM::Matrix3333::operator*(const EigenMatrix3& multi)
{
    Matrix3333 res;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            res.mat[i][j].setZero();
            for (unsigned int k = 0; k < 3; k++) {
                res.mat[i][j] += mat[i][k] * multi(k, j);
            }
        }
    }
    return res;
}

EFEM::Matrix3333 EFEM::Matrix3333::operator*(ScalarType multi)
{
    Matrix3333 res;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            res.mat[i][j] = mat[i][j] * multi;
        }
    }
    return res;
}

EFEM::Matrix3333 EFEM::Matrix3333::transpose()
{
    Matrix3333 res;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            res(i, j) = mat[j][i];
        }
    }
    return res;
}

EFEM::EigenMatrix3 EFEM::Matrix3333::Contract(const EigenMatrix3& multi)
{
    EigenMatrix3 res;
    res.setZero();
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            res += mat[i][j] * multi(i, j);
        }
    }
    return res;
}

EFEM::Matrix3333 EFEM::Matrix3333::Contract(const Matrix3333& multi)
{
    Matrix3333 res;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            res(i, j) = this->Contract(multi.mat[i][j]);
        }
    }
    return res;
}

EFEM::Matrix3333 EFEM::operator*(const EigenMatrix3& multi1, const Matrix3333& multi2)
{
    Matrix3333 res;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            res(i, j).setZero();
            for (unsigned int k = 0; k < 3; k++) {
                res(i, j) += multi1(i, k) * multi2.mat[k][j];
            }
        }
    }
    return res;
}

EFEM::Matrix3333 EFEM::operator*(ScalarType multi1, const Matrix3333& multi2)
{
    Matrix3333 res;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            res(i, j) = multi1 * multi2.mat[i][j];
        }
    }
    return res;
}