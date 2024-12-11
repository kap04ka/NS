#ifndef IEQUATIONSOLVER_H
#define IEQUATIONSOLVER_H

#pragma once
#include <vector>
#include <memory>
#include <string>
#include "ICalculation.h"
#include "ILogger.h"

class IEquationModel
{
public:
    virtual ~IEquationModel() {};
    virtual void solve(std::unique_ptr<ICalculation> calculator, std::shared_ptr<ILogger> logger) = 0;
    virtual void set_params(double u_max, double density, double kinematic_viscosity) = 0;
    virtual void set_time(double tmax, double tau) = 0;
    virtual void set_dims(int nx, int ny, double step) = 0;
    virtual void print_res(const std::string& pathname, std::vector<std::vector<double>> grid) = 0;
    virtual void init_boundary_conditions() = 0;
    virtual void update_boundary_conditions() = 0;
    virtual std::vector<std::vector<double>> get_u() = 0;
    virtual std::vector<std::vector<double>> get_v() = 0;

protected:
    std::shared_ptr<ILogger> logger;
    double u_max;
    double density;
    double kinematic_viscosity;
    double tmax;
    double tau;
    double step;
    int nx;
    int ny;
    std::vector<std::vector<double>> u; // Скорость по оси OX
    std::vector<std::vector<double>> v; // Скорость по оси OY
};

#endif