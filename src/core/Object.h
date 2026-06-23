#ifndef EFEM_CORE_OBJECT_H
#define EFEM_CORE_OBJECT_H

#include "core/IdGenerator.h"
#include <memory>
#include <mutex>

namespace EFEM {
class Object {
public:
    Object()
        : m_id(IdGenerator::instance()->next())
    {
    }

    void assignId(EFEMId id) { m_id = id; }
    EFEMId id() const { return m_id; }

    bool operator==(Object& obj) { return this->m_id == obj.id(); }
    virtual ~Object() = default;

    static std::mutex object_mutex;

private:
    EFEMId m_id;
};
} // namespace EFEM

#endif