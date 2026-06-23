#include "material/Material.h"
#include "math/Math.h"

Eigen::Matrix3d EFEM::NeoHookean::getFirstPiola(Eigen::Matrix3d& F)
{
    double J = F.determinant();
    Eigen::Matrix3d C = F.transpose() * F;
    double I1 = C.trace();
    double I3 = J * J;
    double I1bar = I1 / std::pow(I3, 1.0 / 3.0);
    Eigen::Matrix3d dI1dF = 2 * F;
    Eigen::Matrix3d dJdF = J * F.inverse().transpose();

    // distort
    Eigen::Matrix3d P_I1bar = -2 * I1bar * dJdF / (3.0 * J) + dI1dF / std::pow(I3, 1.0 / 3.0);
    Eigen::Matrix3d P = C10 * P_I1bar;
    // volume
    P += 2 / D1 * (J - 1) * dJdF;
    return P;
}

void EFEM::NeoHookean::initFromJson(const Json& json)
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
    if (json.find("D1") != json.end() && json["D1"].is_number()) {
        this->D1 = json["D1"].template get<double>();
    }
    // t=1.0
    // this->elasticity_coeff = 2.66667 * C10 + 2.0 / D1;
    // t=1.5
    // this->elasticity_coeff = 2.76263 * C10 + 2.0 / D1;
    // t=2.0
    this->elasticity_coeff = 2.94001 * C10 + 2.0 / D1;
    // t=2.5
    // this->elasticity_coeff = 3.13277 * C10 + 2.0 / D1;
    // t=3.0
    // this->elasticity_coeff = 3.32275 * C10 + 2.0 / D1;
    this->type = MaterialType::IN_NEOHOOKEAN;
}