#ifndef ILOGGER_H
#define ILOGGER_H

#include <string>

enum class LogLevel {
    INFO,
    WARNING,
    ERROR
};

class ILogger {
public:
    virtual ~ILogger() = default;
    virtual void log(const std::string& message, LogLevel level) = 0;
    virtual void setLogLevel(LogLevel level) = 0;
};

#endif // ILOGGER_H