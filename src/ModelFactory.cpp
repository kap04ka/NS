#include "ModelFactory.h"
#include <stdexcept>

std::unique_ptr<IEquationModel> ModelFactory::createModel(const std::string& method) {
    if (method == "vorticity-stream") {
        return std::make_unique<VorticityStreamModel>();
    } else if (method == "velocity-pressure") {
        return std::make_unique<VelocityPressureModel>(); 
    } else if (method == "q" || method == "quit") {
        return nullptr;
    } else {
        throw std::invalid_argument("Invalid method.");
    }
}