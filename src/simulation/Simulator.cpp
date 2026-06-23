#include "simulation/Simulator.h"
#include "io/Logger.h"
#include "math/MathType.h"
#include <time.h>

void EFEM::Simulator::initial(std::string& config_file_path)
{
    Eigen::initParallel();
    m_config = std::make_shared<Config>(config_file_path);
    m_model = std::make_shared<Model>(m_config);
    m_constraint = std::make_shared<Constraint>(m_config, m_model);
    if (m_config->outputConfig()->show_log) {
        Logger::info("Simulator initialed");
    }
    if (m_config->outputConfig()->save_frame) {
        m_model->saveVTK("undeformed");
    }
}

void EFEM::Simulator::start()
{
    m_integrator = std::make_shared<Integrator>(m_config, m_model, m_constraint);
    if (m_config->outputConfig()->show_log) {
        Logger::info("Simulator started");
    }

    clock_t stepStartTime, stepEndTime;
    stepStartTime = clock();
    m_integrator->integral();
    m_constraint->updateDx();
    stepEndTime = clock();
    m_config->simulationInfo()->simulate_time += (double)(stepEndTime - stepStartTime) / CLOCKS_PER_SEC;

    if (m_config->outputConfig()->show_log) {
        Logger::info("Simulator step end, use {}s", (double)(stepEndTime - stepStartTime) / CLOCKS_PER_SEC);
    }

    if (m_config->outputConfig()->save_frame) {
        m_model->saveVTK("deformed");
    }
}