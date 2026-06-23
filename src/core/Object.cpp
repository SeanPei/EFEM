#include "core/Object.h"

std::mutex EFEM::Object::object_mutex;
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2