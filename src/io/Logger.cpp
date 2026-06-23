#include "io/Logger.h"

std::mutex EFEM::Logger::mutex;
std::string EFEM::Logger::LOGGER_DISK_LOG_NAME = "EFEM_disk_logger";
std::shared_ptr<EFEM::Logger> EFEM::Logger::pSingleInstance;
std::shared_ptr<EFEM::Logger> EFEM::Logger::get()
{
    if (pSingleInstance.get() == NULL) {
        mutex.lock();
        if (pSingleInstance.get() == NULL) {
            pSingleInstance.reset(new Logger());
        }
        mutex.unlock();
    }
    return pSingleInstance;
}
