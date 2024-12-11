#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#pragma once

#include "ICalculation.h"
#include "ParallelCalculation.h"
#include "SequentialCalculation.h"

#include <memory>
#include <string>

class SolverFactory {
public:
    static std::unique_ptr<ICalculation> createSolver(const std::string& calculator);
};

#endif