#ifndef EFEM_IO_IO_H
#define EFEM_IO_IO_H

#include "Eigen/Dense"
#include "nlohmann/json.hpp"
#include <memory>
#include <string>

using Json = nlohmann::json;

namespace EFEM {
class IO {
public:
    IO(const std::string& config_file_path, const std::string& time_stamp);

    bool checkFile(const std::string& file_path);
    bool checkDir(const std::string& dir_path);
    void ensureDir(const std::string& dir_path);

    std::shared_ptr<Json> configJson() { return m_json; }
    void readFile(std::string file_name,
        const std::function<void(int index, const std::string&)> on_line_callback);
    void loadMesh(std::string file_name,
        const std::function<void(int index, const Eigen::Vector3d&)> on_node_callback,
        const std::function<void(int index, const std::vector<int>&)> on_element_callback);

    std::string workingDir() { return m_working_dir; }
    std::string configFile() { return m_config_file; }
    std::string logFile() { return m_log_file; }

private:
    std::shared_ptr<Json> m_json;
    std::string m_working_dir = "";
    std::string m_config_file = "";
    std::string m_log_file = "";
};
}

#endif