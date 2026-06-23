#ifndef EFEM_MATERIAL_MATERIAL_H
#define EFEM_MATERIAL_MATERIAL_H

#include "math/Math.h"
#include "nlohmann/json.hpp"

using Json = nlohmann::json;
namespace EFEM {
enum class MaterialType {
    NEOHOOKEAN,
    MOONEY_RIVLIN,
    YEOH,
    IN_NEOHOOKEAN,
    IN_MOONEY_RIVLIN,
    IN_YEOH,
};

class Material {
#define J_CLAMP_VALUE 0.1
public:
    int id = -1;
    MaterialType type;
    double density = 0.f;
    double elasticity_coeff = 0.f;

    virtual Eigen::Matrix3d getFirstPiola(Eigen::Matrix3d& F) = 0;
    virtual void initFromJson(const Json& json) = 0;
};

class NeoHookean : public Material {
public:
    double C10 = 0, D1 = 0;

    Eigen::Matrix3d getFirstPiola(Eigen::Matrix3d& F);
    void initFromJson(const Json& json);
};

class MooneyRivlin : public Material {
public:
    double C10 = 0, C01 = 0, D1 = 0;

    Eigen::Matrix3d getFirstPiola(Eigen::Matrix3d& F);
    void initFromJson(const Json& json);
};

class Yeoh : public Material {
public:
    double C10 = 0, C20 = 0, C30 = 0, D1 = 0;

    Eigen::Matrix3d getFirstPiola(Eigen::Matrix3d& F);
    void initFromJson(const Json& json);
};

}

#endif