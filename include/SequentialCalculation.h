#ifndef SEQUENTIALCALCULATION_H
#define SEQUENTIALCALCULATION_H

#pragma once
#include "ICalculation.h"

class SequentialCalculation : public ICalculation { 
public:
    SequentialCalculation() = default;
    ~SequentialCalculation() override = default;

    void jacobi(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& extends, double epsilon, double step) override;
    void four_point_pattern(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid, double step) override;
    void reverse_four_point_pattern(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid, double step) override;
    void fill_zeros_vector_1D(std::vector<double>& grid) override;
    void fill_zeros_vector_2D(std::vector<std::vector<double>>& grid) override;
    void copy_vector(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid) override;
    void div_calculation(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid, const std::vector<std::vector<double>>& temp_grid_second, double step) override;
    void helmholtz_solve(std::vector<std::vector<double>>& grid, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, double step, double tau, double kinematic_viscosity) override;
    void update_velocity(std::vector<std::vector<double>>& grid, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& pressure, double step, double tau, double kinematic_viscosity, double density) override;
    void update_pressure(std::vector<std::vector<double>> &grid, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, double step, double c, double tau) override;
};

#endif