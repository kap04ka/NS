#ifndef VORTICITYSTEAMMODEL_H
#define VORTICITYSTEAMMODEL_H

#pragma once
#include "IEquationModel.h"

class VorticityStreamModel : public IEquationModel
{
public:
    VorticityStreamModel();
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
    double Q;
    double epsilon;
    std::vector<std::vector<double>> vorticity;
    std::vector<std::vector<double>> vorticity_new;
    std::vector<std::vector<double>> stream_function;
    std::vector<std::vector<double>> stream_function_new;
    std::vector<std::vector<double>> obstruction;
};
#endif