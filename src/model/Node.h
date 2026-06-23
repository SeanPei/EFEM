#ifndef EFEM_MODEL_NODE_H
#define EFEM_MODEL_NODE_H
#include "Eigen/Dense"
#include "core/Object.h"
#include <memory>
#include <unordered_map>
#include <vector>

namespace EFEM {
class Element;

class Node : public Object {
public:
    int index;
    std::vector<std::shared_ptr<Element>> elements;
    Eigen::Vector3d position;
    Eigen::Vector3d dx;
    std::unordered_map<std::string, std::shared_ptr<Object>> extra_info;
};
}

#endif
