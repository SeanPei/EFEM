#ifndef EFEM_IO_CONFIG_H
#define EFEM_IO_CONFIG_H
#include "Eigen/Dense"
#include "io/IO.h"
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace EFEM {

struct SimulationInfo {
    std::string time_stamp = "";
    double simulate_time = 0;
};

struct AnalysisConfig {
    // 求解器算法设置项
    enum NLSolver {
        LBFGS_AH,
        LBFGS_I
    };
    NLSolver nl_solver = NLSolver::LBFGS_AH;
    int max_solver_iters = 0; // 最大迭代步数
    int max_line_search = 20; // 最大线性搜索步数
    double dx_criterion = 1e-8; // dx终止判据
    double residual_criterion = 1e-6; // residual终止判据
    // 全局设置
    enum Gravity {
        NON,
        X,
        NX,
        Y,
        NY,
        Z,
        NZ
    };
    Gravity gravity = Gravity::NON;
    bool consistent_inertia = false; // 一致质量矩阵
};

struct OutputConfig {
    bool save_frame = true;
    // 日志选项
    bool disk_log = true; // 是否将log写入磁盘
    bool show_log = true;
};

class Config {
public:
    Config(const std::string& config_file_path);
    std::shared_ptr<IO> io() { return m_io; }
    std::string workingDir();
    std::string configFile();
    std::string logFile();

    std::shared_ptr<SimulationInfo> simulationInfo() { return m_simulation_info; }
    std::shared_ptr<AnalysisConfig> analysisConfig() { return m_analysis_config; }
    std::shared_ptr<OutputConfig> outputConfig() { return m_output_config; }

private:
    std::shared_ptr<SimulationInfo> m_simulation_info;
    std::shared_ptr<IO> m_io;
    std::shared_ptr<AnalysisConfig> m_analysis_config;
    std::shared_ptr<OutputConfig> m_output_config;
};
}

#endif