#include "exceptions.hpp"

ImageException::ImageException(const std::string& message)
    : std::runtime_error(message) {}
