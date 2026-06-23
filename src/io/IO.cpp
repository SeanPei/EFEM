#include "io/IO.h"
#include "io/Logger.h"
#include "io/MeshParser.h"
#include "spdlog/spdlog.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <fstream>

bool EFEM::IO::checkFile(const std::string& file_path)
{
    boost::filesystem::path file(file_path);
    return boost::filesystem::exists(file) && boost::filesystem::is_regular_file(file);
}

bool EFEM::IO::checkDir(const std::string& dir_path)
{
    boost::filesystem::path dir(dir_path);
    return boost::filesystem::exists(dir) && boost::filesystem::is_directory(dir);
}

void EFEM::IO::ensureDir(const std::string& dir_path)
{
    if (!checkDir(dir_path)) {
        boost::filesystem::create_directories(dir_path);
    }
}

EFEM::IO::IO(const std::string& config_file_path, const std::string& time_stamp)
{
    std::string path_separator;
#ifdef WIN32
    path_separator = "/";
#else
    path_separator = "/";
#endif
    m_working_dir = config_file_path.substr(0, config_file_path.find_last_of(path_separator) + 1);
    m_config_file = config_file_path.substr(config_file_path.find_last_of(path_separator) + 1);

    // 打开配置文件
    std::string file = m_working_dir + m_config_file;
    if (!checkFile(file)) {
        spdlog::error("config file {} not exist!", file);
        return;
    }
    std::ifstream ifs(file);
    try {
        if (!ifs.is_open())
            throw std::runtime_error(
                std::string("Input file not found in specified location: ") + file);
    } catch (const std::runtime_error& except) {
        spdlog::error("{}", except.what());
        std::terminate();
    }

    try {
        m_json = std::make_shared<Json>(Json::parse(ifs));
    } catch (const std::runtime_error& except) {
        spdlog::error("{}", except.what());
        std::terminate();
    }

    std::string log_dir = m_working_dir + "log/";
    ensureDir(log_dir);
    m_log_file = log_dir + time_stamp + ".txt";
}

void EFEM::IO::readFile(std::string file_name,
    const std::function<void(int index, const std::string&)> on_line_callback)
{
    std::string file_path = m_working_dir + file_name;
    if (!checkFile(file_path)) {
        Logger::error("data file {} not exist!", file_name);
        return;
    }
    std::fstream file;
    file.open(file_path.c_str(), std::ios::in);
    try {
        if (file.is_open() && file.good()) {
            // Line
            std::string line;
            int index = 0;
            while (std::getline(file, line)) {
                boost::algorithm::trim(line);
                // ignore comment lines (# or !) or blank lines
                if ((line.find('#') == std::string::npos) && (line.find('!') == std::string::npos) && (line != "")) {
                    on_line_callback(index, line);
                    index++;
                }
            }
        }
        file.close();
    } catch (std::exception& exception) {
        Logger::error("Read element stresses: {}", exception.what());
        file.close();
    }
}

void EFEM::IO::loadMesh(std::string file_name,
    const std::function<void(int index, const Eigen::Vector3d&)> on_node_callback,
    const std::function<void(int index, const std::vector<int>&)> on_element_callback)
{
    std::string file_path = m_working_dir + file_name;
    MeshParser meshParser;
    meshParser.parse(file_path, on_node_callback, on_element_callback);
}