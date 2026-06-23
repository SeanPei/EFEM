#include "material/Material.h"
#include "math/Math.h"

Eigen::Matrix3d EFEM::MooneyRivlin::getFirstPiola(Eigen::Matrix3d& F)
{
    Eigen::Matrix3d C = F.transpose() * F;
    double J = F.determinant();
    double I1 = C.trace();
    double I2 = 0.5 * (std::pow(I1, 2) - (C.transpose() * C).trace());
    double I3 = J * J;
    double I1bar = I1 / std::pow(I3, 1.0 / 3.0);
    double I2bar = I2 / std::pow(I3 * I3, 1.0 / 3.0);
    Eigen::Matrix3d dI1dF = 2 * F;
    Eigen::Matrix3d dI2dF = 2 * F * (I1 * Eigen::Matrix3d::Identity() - C);
    Eigen::Matrix3d dJdF = J * F.inverse().transpose();

    // distort
    Eigen::Matrix3d P_I1bar = -2 * I1bar * dJdF / (3.0 * J) + dI1dF / std::pow(I3, 1.0 / 3.0);
    Eigen::Matrix3d P_I2bar = -4 * I2bar * dJdF / (3.0 * J) + dI2dF / std::pow(I3 * I3, 1.0 / 3.0);
    Eigen::Matrix3d P = C10 * P_I1bar + C01 * P_I2bar;
    // volume
    P += 2 / D1 * (J - 1) * dJdF;
    return P;
}

void EFEM::MooneyRivlin::initFromJson(const Json& json)
{
    if (json.find("id") != json.end() && json["id"].is_number()) {
        this->id = json["id"].template get<int>();
    }
    if (json.find("density") != json.end() && json["density"].is_number()) {
        this->density = json["density"].template get<double>();
    }
    if (json.find("C10") != json.end() && json["C10"].is_number()) {
        this->C10 = json["C10"].template get<double>();
    }
    if (json.find("C01") != json.end() && json["C01"].is_number()) {
        this->C01 = json["C01"].template get<double>();
    }
    if (json.find("D1") != json.end() && json["D1"].is_number()) {
        this->D1 = json["D1"].template get<double>();
    }
    // t=1.0
    // this->elasticity_coeff = 2.66667 * C01 + 2.66667 * C10 + 2.0 / D1;
    // t=1.5
    // this->elasticity_coeff = 3.0575 * C01 + 2.76263 * C10 + 2.0 / D1;
    // t=2.0
    this->elasticity_coeff = 3.818 * C01 + 2.94001 * C10 + 2.0 / D1;
    // t=2.5
    // this->elasticity_coeff = 4.70324 * C01 + 3.13277 * C10 + 2.0 / D1;
    // t=3.0
    // this->elasticity_coeff = 5.63899 * C01 + 3.32275 * C10 + 2.0 / D1;
    this->type = MaterialType::IN_MOONEY_RIVLIN;
}