#ifndef EFEM_SIMULATION_SIMULATOR_H
#define EFEM_SIMULATION_SIMULATOR_H

#include "fe/Constraint.h"
#include "io/Config.h"
#include "model/Model.h"
#include "simulation/Integrator.h"
#include <memory>

namespace EFEM {
class Simulator {
public:
    void initial(std::string& config_file_path);
    void start();
    std::string info();

    std::shared_ptr<Config> getConfig() { return m_config; }
    std::shared_ptr<Model> getModel() { return m_model; }
    std::shared_ptr<Constraint> getConstraint() { return m_constraint; }

private:
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
    std::shared_ptr<Constraint> m_constraint;
    std::shared_ptr<Integrator> m_integrator;
};
}

#endif