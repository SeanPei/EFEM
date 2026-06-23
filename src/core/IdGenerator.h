#ifndef EFEM_CORE_IDGENERATOR_H
#define EFEM_CORE_IDGENERATOR_H

namespace EFEM {

using EFEMId = unsigned long long;

class IdGenerator {
public:
    static IdGenerator* instance()
    {
        static IdGenerator idGenerator;
        return &idGenerator;
    }

    EFEMId next() { return _id++; }

private:
    IdGenerator()
        : _id(0)
    {
    }

    EFEMId _id;
};
}
#endif