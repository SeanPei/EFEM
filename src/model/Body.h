#ifndef EFEM_MODEL_BODY_H
#define EFEM_MODEL_BODY_H

#include "io/Config.h"
#include "model/Element.h"
#include "model/Node.h"
#include <memory>
#include <vector>

namespace EFEM {

enum BodyType {
    SOFT,
    RIGID
};

class Body {
public:
    int id = -1;
    Body(const int id, const std::shared_ptr<Config> config, const std::vector<std::pair<std::string, std::vector<double>>>& poses, const std::string& mesh_file);
    double volume();

    std::vector<std::shared_ptr<Element>>& getElements() { return m_elements; }
    std::vector<std::shared_ptr<Node>>& getNodes() { return m_nodes; }

private:
    std::vector<std::shared_ptr<Element>> m_elements;
    std::vector<std::shared_ptr<Node>> m_nodes;
    std::vector<std::shared_ptr<Node>> m_render_nodes;

    void buildNode(int index, const Eigen::Vector3d& node_pos);
    void buildElement(int index, const std::vector<int>& node_list);
};
}

#endif