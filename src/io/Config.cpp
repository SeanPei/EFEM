#include "io/Config.h"
#include "io/Logger.h"
#include "material/Material.h"
#include <ctime>

EFEM::Config::Config(const std::string& config_file_path)
{
    // 初始化simulation参数
    m_simulation_info = std::make_shared<EFEM::SimulationInfo>();
    std::time_t t = std::time(0); // get time now
    std::tm* now = std::localtime(&t);
    m_simulation_info->time_stamp = std::to_string(now->tm_year + 1900) + "-" + std::to_string(now->tm_mon + 1) + "-" + std::to_string(now->tm_mday) + "-" + std::to_string(now->tm_hour) + "-" + std::to_string(now->tm_min) + "-" + std::to_string(now->tm_sec);

    // 初始化io
    m_io = std::make_shared<EFEM::IO>(config_file_path, m_simulation_info->time_stamp);

    // 初始化logger
    Logger::initial(this);

    auto config = m_io->configJson();
    // 解析analysis参数
    m_analysis_config = std::make_shared<EFEM::AnalysisConfig>();
    if (config->find("analysis") != config->end()) {
        Json analysis = *(config->find("analysis"));
        if (analysis.find("nl_solver") != analysis.end() && analysis["nl_solver"].is_string()) {
            std::string nl_solver = analysis["nl_solver"].template get<std::string>();
            if (nl_solver == "lbfgs_ah") {
                m_analysis_config->nl_solver = AnalysisConfig::NLSolver::LBFGS_AH;
            } else if (nl_solver == "lbfgs") {
                m_analysis_config->nl_solver = AnalysisConfig::NLSolver::LBFGS_I;
            }
        }
        if (analysis.find("max_solver_iters") != analysis.end() && analysis["max_solver_iters"].is_number()) {
            m_analysis_config->max_solver_iters = analysis["max_solver_iters"].template get<int>();
        }
        if (analysis.find("max_line_search") != analysis.end() && analysis["max_line_search"].is_number()) {
            m_analysis_config->max_line_search = analysis["max_line_search"].template get<int>();
        }
        if (analysis.find("dx_criterion") != analysis.end() && analysis["dx_criterion"].is_number()) {
            m_analysis_config->dx_criterion = analysis["dx_criterion"].template get<double>();
        }
        if (analysis.find("residual_criterion") != analysis.end() && analysis["residual_criterion"].is_number()) {
            m_analysis_config->residual_criterion = analysis["residual_criterion"].template get<double>();
        }
        if (analysis.find("gravity") != analysis.end() && analysis["gravity"].is_string()) {
            std::string gravity = analysis["gravity"].template get<std::string>();
            if (gravity == "non") {
                m_analysis_config->gravity = AnalysisConfig::Gravity::NON;
            } else if (gravity == "x") {
                m_analysis_config->gravity = AnalysisConfig::Gravity::X;
            } else if (gravity == "-x") {
                m_analysis_config->gravity = AnalysisConfig::Gravity::NX;
            } else if (gravity == "y") {
                m_analysis_config->gravity = AnalysisConfig::Gravity::Y;
            } else if (gravity == "-y") {
                m_analysis_config->gravity = AnalysisConfig::Gravity::NY;
            } else if (gravity == "z") {
                m_analysis_config->gravity = AnalysisConfig::Gravity::Z;
            } else if (gravity == "-z") {
                m_analysis_config->gravity = AnalysisConfig::Gravity::NZ;
            }
        }
        if (analysis.find("consistent_inertia") != analysis.end() && analysis["consistent_inertia"].is_boolean()) {
            m_analysis_config->consistent_inertia = analysis["consistent_inertia"].template get<bool>();
        }
    }

    // 解析output参数
    m_output_config = std::make_shared<EFEM::OutputConfig>();
    if (config->find("output") != config->end()) {
        Json output = *(config->find("output"));
        if (output.find("save_frame") != output.end() && output["save_frame"].is_boolean()) {
            m_output_config->save_frame = output["save_frame"].template get<bool>();
        }
        if (output.find("show_log") != output.end() && output["show_log"].is_boolean()) {
            m_output_config->show_log = output["show_log"].template get<bool>();
        }
        if (output.find("disk_log") != output.end() && output["disk_log"].is_boolean()) {
            m_output_config->disk_log = output["disk_log"].template get<bool>();
        }
    }
}

std::string EFEM::Config::workingDir() { return m_io ? m_io->workingDir() : ""; }
std::string EFEM::Config::configFile() { return m_io ? m_io->configFile() : ""; }
std::string EFEM::Config::logFile() { return m_io ? m_io->logFile() : ""; }
