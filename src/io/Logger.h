#ifndef EFEM_IO_LOGGER_H
#define EFEM_IO_LOGGER_H

#include "io/Config.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/spdlog.h"
#include <mutex>

namespace EFEM {
class Logger {
public:
    static inline void initial(Config* config)
    {
        Logger::get()->initialImpl(config);
    }

    template <typename... Args>
    static inline void info(std::string fmt, const Args&... args)
    {
        Logger::get()->infoImpl(fmt, args...);
    }

    template <typename... Args>
    static inline void warn(std::string fmt, const Args&... args)
    {
        Logger::get()->warnImpl(fmt, args...);
    }

    template <typename... Args>
    static inline void error(std::string fmt, const Args&... args)
    {
        Logger::get()->errorImpl(fmt, args...);
    }

private:
    Logger() = default;
    static std::shared_ptr<Logger> get();
    static std::mutex mutex;
    static std::shared_ptr<Logger> pSingleInstance;
    static std::string LOGGER_DISK_LOG_NAME;

    Config* m_config_ptr;

    void initialImpl(Config* config)
    {
        m_config_ptr = config;
    }

    template <typename... Args>
    void infoImpl(std::string fmt, const Args&... args)
    {
        if (m_config_ptr->outputConfig()->disk_log && initial_disk_log()) {
            spdlog::get(LOGGER_DISK_LOG_NAME)->info(fmt, args...);
            spdlog::get(LOGGER_DISK_LOG_NAME)->flush();
        }
        if (m_config_ptr->outputConfig()->show_log) {
            spdlog::info(fmt, args...);
        }
    }

    template <typename... Args>
    void warnImpl(std::string fmt, const Args&... args)
    {
        if (m_config_ptr->outputConfig()->disk_log && initial_disk_log()) {
            spdlog::get(LOGGER_DISK_LOG_NAME)->warn(fmt, args...);
            spdlog::get(LOGGER_DISK_LOG_NAME)->flush();
        }
        if (m_config_ptr->outputConfig()->show_log) {
            spdlog::warn(fmt, args...);
        }
    }

    template <typename... Args>
    void errorImpl(std::string msg, const Args&... args)
    {
        if (m_config_ptr->outputConfig()->disk_log && initial_disk_log()) {
            spdlog::get(LOGGER_DISK_LOG_NAME)->error(msg, args...);
            spdlog::get(LOGGER_DISK_LOG_NAME)->flush();
        }
        if (m_config_ptr->outputConfig()->show_log) {
            spdlog::error(msg, args...);
        }
    }

    bool initial_disk_log()
    {
        static bool inited = false;
        static std::string last_log_file = "";
        std::string log_file = m_config_ptr->logFile();
        if (log_file.length() == 0) {
            return false;
        }
        inited = inited & last_log_file == log_file;
        if (!inited) {
            mutex.lock();
            if (!inited) {
                inited = true;
                if (spdlog::get(LOGGER_DISK_LOG_NAME) != nullptr) {
                    spdlog::drop(LOGGER_DISK_LOG_NAME);
                }
                spdlog::basic_logger_mt(LOGGER_DISK_LOG_NAME, log_file);
                last_log_file = log_file;
            }
            mutex.unlock();
        }
        return true;
    }
};
}

#endif