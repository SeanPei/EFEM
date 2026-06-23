#ifndef EFEM_IO_MESH_PARSER_H
#define EFEM_IO_MESH_PARSER_H

#include "Eigen/Dense"
#include "io/Config.h"
#include "io/Logger.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

namespace EFEM {

class MeshParser {
public:
    MeshParser()
    {
    }

    void parse(std::string file_path,
        const std::function<void(int index, const Eigen::Vector3d&)> on_node_callback,
        const std::function<void(int index, const std::vector<int>&)> on_element_callback)
    {
        if (!checkFile(file_path)) {
            Logger::error("data file {} not exist!", file_path);
            return;
        }
        int nnode = 0;
        int nelement = 0;
        std::fstream file;
        file.open(file_path.c_str(), std::ios::in);
        try {
            if (file.is_open() && file.good()) {
                // Line
                std::string line;
                bool available = nextLine(file, line);
                while (available) {
                    if (std::strncmp(line.c_str(), "#node", 5) == 0) {
                        available = nextLine(file, line);
                        while (available && !(std::strncmp(line.c_str(), "#", 1) == 0)) {
                            std::istringstream istream(line);
                            Eigen::Vector3d node;
                            for (int i = 0; i < 3; ++i) {
                                istream >> node(i);
                            }
                            on_node_callback(nnode, node);
                            nnode++;
                            available = nextLine(file, line);
                        }
                        continue;
                    }

                    if (std::strncmp(line.c_str(), "#element", 8) == 0) {
                        available = nextLine(file, line);
                        while (available && !(std::strncmp(line.c_str(), "#", 1) == 0)) {
                            std::vector<std::string> s_node_list;
                            boost::split(s_node_list, line, boost::is_any_of(", "), boost::token_compress_on);
                            std::vector<int> node_list;
                            for (auto s_node : s_node_list) {
                                node_list.push_back(std::stoi(s_node));
                            }
                            on_element_callback(nelement, node_list);
                            nelement++;
                            available = nextLine(file, line);
                        }
                        continue;
                    }

                    available = nextLine(file, line);
                }
            }
            file.close();
        } catch (std::exception& exception) {
            Logger::error("Read element stresses: {}", exception.what());
            file.close();
        }
    }

    bool nextLine(std::fstream& file, std::string& line)
    {
        bool available = !file.eof();
        if (available) {
            std::getline(file, line);
            boost::algorithm::trim(line);
            while (std::strncmp(line.c_str(), "##", 2) == 0 || std::strncmp(line.c_str(), "!!", 2) == 0 || (line == "")) {
                available = !file.eof();
                if (!available) {
                    break;
                }
                std::getline(file, line);
                boost::algorithm::trim(line);
            }
        }
        if (!available) {
            line = "";
        }
        return available;
    }

    bool checkFile(const std::string& file_path)
    {
        boost::filesystem::path file(file_path);
        return boost::filesystem::exists(file) && boost::filesystem::is_regular_file(file);
    }
};
}

#endif