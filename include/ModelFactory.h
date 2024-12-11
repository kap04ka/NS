#ifndef MODELFACTORY_H
#define MODELFACTORY_H

#pragma once

#include "IEquationModel.h"
#include "VelocityPressureModel.h"
#include "VorticityStreamModel.h"

#include <memory>
#include <string>

class ModelFactory {
public:
    static std::unique_ptr<IEquationModel> createModel(const std::string& method);
};

#endif