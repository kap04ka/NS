#ifndef VELOCITYPRESSUREMODEL_H
#define VELOCITYPRESSUREMODEL_H

#pragma once
#include "IEquationModel.h"

class VelocityPressureModel : public IEquationModel
{
public:
    VelocityPressureModel();
    void solve(std::unique_ptr<ICalculation> calculator, std::shared_ptr<ILogger> logger) override;
    void set_params(double u_max, double density, double kinematic_viscosity) override;
    void set_time(double tmax, double tau) override;
    void set_dims(int nx, int ny, double step) override;
    void print_res(const std::string& pathname, std::vector<std::vector<double>> grid) override;
    void init_boundary_conditions() override;
    void update_boundary_conditions() override;

    std::vector<std::vector<double>> get_u() {
        return u;
    }

    std::vector<std::vector<double>> get_v() {
        return v;
    }

private:
    double c;
    std::vector<std::vector<double>> pressure;
    std::vector<std::vector<double>> div_velocity;
    std::vector<std::vector<double>> u_new;
    std::vector<std::vector<double>> v_new;
};

#endif