#ifndef FLUIDDYNAMIC_H
#define FLUIDDYNAMIC_H

#pragma once
#include "IEquationModel.h"
#include "ILogger.h"
#include <string>

class FlyidDynamic
{
public:
    FlyidDynamic(std::unique_ptr<IEquationModel> model, std::unique_ptr<ICalculation> calculator, std::shared_ptr<ILogger> logger);

    void set_dim_params(double height, double width, double step);
    void set_time_params(double tmax, double tau);
    void set_model_params(double u_max, double density, double kinematic_viscosity);
    void start();
    void print_results(const std::string& pathname_u, const std::string& pathname_v);

protected:
    std::shared_ptr<ILogger> logger;
    std::unique_ptr<IEquationModel> model;
    std::unique_ptr<ICalculation> calculator;
    double nx;
    double ny;
    double step;
    double tmax;
    double tau;
    double u_max;
    double density;
    double kinematic_viscosity;
};

#endif