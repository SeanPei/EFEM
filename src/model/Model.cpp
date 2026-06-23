#include "model/Model.h"
#include <chrono>
#include <ctime>
#include <vtkPointData.h>
#include <vtkUnstructuredGridWriter.h>

EFEM::Model::Model(std::shared_ptr<Config> config)
{
    m_config = config;
    std::shared_ptr<Json> config_json = config->io()->configJson();
    if (!(config_json->find("model") != config_json->end() && config_json->find("model")->is_array())) {
        return;
    }
    std::vector<Json> model_jsons = config_json->find("model")->template get<std::vector<Json>>();
    for (auto mod_json : model_jsons) {
        if (mod_json["type"].is_string()) {
            std::vector<std::pair<std::string, std::vector<double>>> poses;
            if (mod_json.find("pose") != mod_json.end() && mod_json["pose"].is_array()) {
                std::vector<Json> mod_pose_jsons = mod_json.find("pose")->template get<std::vector<Json>>();
                for (auto mod_pose_json : mod_pose_jsons) {
                    for (auto mod_pos : mod_pose_json.items()) {
                        poses.push_back({ mod_pos.key(), mod_pos.value().template get<std::vector<double>>() });
                    }
                }
            }
            int id = -1;
            if (mod_json.find("id") != mod_json.end() && mod_json["id"].is_number()) {
                id = mod_json["id"].template get<int>();
            }
            std::string mesh_file = mod_json["mesh_file"].template get<std::string>();

            if (mod_json.find("mesh_file") != mod_json.end() && mod_json["mesh_file"].is_string()) {
                std::shared_ptr<Body> body = std::make_shared<Body>(id, config, poses, mesh_file);
                m_bodies.push_back(body);
                m_body_map[id] = body;
            }
        }
    }
}

double EFEM::Model::volume()
{
    double volume = 0;
    for (auto body : m_bodies) {
        volume += body->volume();
    }
    return volume;
}

void EFEM::Model::saveVTK(std::string name)
{
    vtkSmartPointer<vtkUnstructuredGrid> scene_data = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> scene_points = vtkSmartPointer<vtkPoints>::New();
    std::vector<std::shared_ptr<Body>> bodies = this->getBodies();
    for (auto body : bodies) {
        std::vector<std::shared_ptr<Node>> nodes = body->getNodes();
        for (int i = 0; i < nodes.size(); i++) {
            Eigen::Vector3d position = nodes[i]->position + nodes[i]->dx;
            scene_points->InsertNextPoint(position.data());
        }
    }
    scene_data->SetPoints(scene_points);

    int body_node_offset = 0;
    for (auto body : bodies) {
        std::vector<std::shared_ptr<Node>> nodes = body->getNodes();
        std::vector<std::shared_ptr<Element>> elements = body->getElements();
        for (int i = 0; i < elements.size(); i++) {
            vtkSmartPointer<vtkIdList> cell_point_lists = vtkSmartPointer<vtkIdList>::New();
            int e_node_nums = elements[i]->node_indexs.size();
            for (int j = 0; j < e_node_nums; j++) {
                cell_point_lists->InsertNextId(elements[i]->node_indexs[j] + body_node_offset);
            }
            VTKCellType cell_type;
            switch (e_node_nums) {
            case 4:
                cell_type = VTKCellType::VTK_TETRA;
                break;
            case 10:
                cell_type = VTKCellType::VTK_QUADRATIC_TETRA;
                break;
            case 8:
                cell_type = VTKCellType::VTK_HEXAHEDRON;
                break;
            case 20:
                cell_type = VTKCellType::VTK_QUADRATIC_HEXAHEDRON;
                break;
            }
            scene_data->InsertNextCell(cell_type, cell_point_lists);
        }
        body_node_offset += nodes.size();
    }

    auto tNow = std::chrono::system_clock::now();
    auto tSeconds = std::chrono::duration_cast<std::chrono::seconds>(tNow.time_since_epoch());
    auto secNow = tSeconds.count();
    tm tmNow;
#ifdef _WIN32
    localtime_s(&tmNow, &secNow);
#elif __linux__
    localtime_r(&secNow, &tmNow);
#endif
    std::ostringstream oss_dir;
    oss_dir << std::put_time(&tmNow, "%Y-%m-%d");
    std::ostringstream oss_file;
    oss_file << std::put_time(&tmNow, "%Y-%m-%d-%H-%M-%S");
    auto tMilli = std::chrono::duration_cast<std::chrono::milliseconds>(tNow.time_since_epoch());
    auto ms = tMilli - tSeconds;
    oss_file << "-" << std::setfill('0') << std::setw(3) << ms.count();
    std::string vtk_dir = m_config->workingDir() + "results/" + oss_dir.str() + "/";
    m_config->io()->ensureDir(vtk_dir);
    std::string path = vtk_dir + oss_file.str() + "-" + name + ".vtk";
    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(path.c_str());
    writer->SetInputData(scene_data);
    writer->Write();
}