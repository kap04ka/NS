#ifndef ICALCULATION_H 
#define ICALCULATION_H

#pragma once
#include <vector>

class ICalculation {
public:
    virtual ~ICalculation() = default;

    virtual void jacobi(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& extends, double epsilon, double step);
    virtual void four_point_pattern(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid, double step);
    virtual void reverse_four_point_pattern(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid, double step);
    virtual void fill_zeros_vector_1D(std::vector<double>& grid);
    virtual void fill_zeros_vector_2D(std::vector<std::vector<double>>& grid);
    virtual void copy_vector(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid);
    virtual void div_calculation(std::vector<std::vector<double>>& grid, const std::vector<std::vector<double>>& temp_grid, const std::vector<std::vector<double>>& temp_grid_second, double step);
    virtual void helmholtz_solve(std::vector<std::vector<double>>& grid, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, double step, double tau, double kinematic_viscosity);
    virtual void update_velocity(std::vector<std::vector<double>>& grid, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& pressure, double step, double tau, double kinematic_viscosity, double density);
    virtual void update_pressure(std::vector<std::vector<double>> &grid, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, double step, double c, double tau);

};

#endif