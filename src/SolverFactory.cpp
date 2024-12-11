#include "SolverFactory.h"
#include <stdexcept>

std::unique_ptr<ICalculation> SolverFactory::createSolver(const std::string& calculator) {
    if (calculator == "parallel") {
        return std::make_unique<ParallelCalculation>();
    } else if (calculator == "sequential") {
        return std::make_unique<SequentialCalculation>(); 
    } else if (calculator == "q" || calculator == "quit") {
        return nullptr;
    } else {
        throw std::invalid_argument("Invalid calculator.");
    }
}