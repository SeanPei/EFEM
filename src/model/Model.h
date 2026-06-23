#ifndef EFEM_MODEL_MODEL_H
#define EFEM_MODEL_MODEL_H

#include "model/Body.h"

namespace EFEM {
class Model {
public:
    Model(std::shared_ptr<Config> config);
    void reset();

    std::shared_ptr<Body> getBody(int id) { return m_body_map[id]; }
    std::vector<std::shared_ptr<Body>> getBodies() { return m_bodies; }
    double volume();
    void saveVTK(std::string name);

private:
    std::shared_ptr<Config> m_config;
    std::vector<std::shared_ptr<Body>> m_bodies;
    std::unordered_map<int, std::shared_ptr<Body>> m_body_map;
};
}

#endif