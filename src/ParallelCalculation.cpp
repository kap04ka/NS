#include "ParallelCalculation.h"
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <iostream> 

void ParallelCalculation::copy_vector(std::vector<std::vector<double>> &grid, const std::vector<std::vector<double>> &temp_grid)
{
#pragma omp parallel for collapse(2)
    for (int i = 1; i < grid.size() - 1; ++i)
    {
        for (int j = 1; j < grid[0].size() - 1; ++j)
        {
            grid[i][j] = temp_grid[i][j];
        }
    }
}

void ParallelCalculation::fill_zeros_vector_1D(std::vector<double> &grid)
{
    std::fill(grid.begin(), grid.end(), 0.0);
}

void ParallelCalculation::fill_zeros_vector_2D(std::vector<std::vector<double>> &grid)
{
    for (auto &row : grid)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
}

void ParallelCalculation::jacobi(std::vector<std::vector<double>> &grid, const std::vector<std::vector<double>> &extends, double epsilon, double step)
{
    double eps{};
    int rows = grid.size();
    int cols = grid[0].size();

    std::vector<std::vector<double>> temp_grid(rows, std::vector<double>(cols, 0.0));
    do
    {
        eps = -1.0;
#pragma omp parallel for reduction(max : eps) collapse(2)
        for (int i = 1; i < rows - 1; ++i)
        {
            for (int j = 1; j < cols - 1; ++j)
            {
                temp_grid[i][j] = (step * step * extends[i][j] + grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1]);
                eps = std::max(eps, fabs(temp_grid[i][j] - grid[i][j]));
            }
        }
        copy_vector(grid, temp_grid);
    } while (eps >= epsilon);
}

void ParallelCalculation::four_point_pattern(std::vector<std::vector<double>> &grid, const std::vector<std::vector<double>> &temp_grid, double step)
{
    int rows = grid.size();
    int cols = grid[0].size();

#pragma omp parallel for collapse(2)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j = 1; j < cols - 1; ++j)
        {
            grid[i][j] = (temp_grid[i + 1][j + 1] + temp_grid[i - 1][j + 1] - temp_grid[i + 1][j - 1] - temp_grid[i - 1][j - 1]) / (4 * step);
        }
    }
}

void ParallelCalculation::reverse_four_point_pattern(std::vector<std::vector<double>> &grid, const std::vector<std::vector<double>> &temp_grid, double step)
{
    int rows = grid.size();
    int cols = grid[0].size();

#pragma omp parallel for collapse(2)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j = 1; j < cols - 1; ++j)
        {
            grid[i][j] = -(temp_grid[i + 1][j + 1] - temp_grid[i - 1][j + 1] + temp_grid[i + 1][j - 1] - temp_grid[i - 1][j - 1]) / (4 * step);
        }
    }
}

void ParallelCalculation::div_calculation(std::vector<std::vector<double>> &grid, const std::vector<std::vector<double>> &temp_grid, const std::vector<std::vector<double>> &temp_grid_second, double step)
{
    int rows = grid.size();
    int cols = grid[0].size();

#pragma omp parallel for collapse(2)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j = 1; j < cols - 1; ++j)
        {
            grid[i][j] = ((temp_grid[i + 1][j + 1] + temp_grid[i + 1][j - 1]) - (temp_grid[i - 1][j - 1] + temp_grid[i - 1][j + 1])) + ((temp_grid_second[i - 1][j + 1] + temp_grid_second[i + 1][j + 1]) - (temp_grid_second[i - 1][j - 1] + temp_grid_second[i + 1][j - 1])) / (4 * step);
        }
    }
}

void ParallelCalculation::helmholtz_solve(std::vector<std::vector<double>> &grid, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, double step, double tau, double kinematic_viscosity)
{
    int rows = grid.size();
    int cols = grid[0].size();

    std::vector<std::vector<double>> temp_grid(rows, std::vector<double>(cols, 0.0));

#pragma omp parallel for collapse(2)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j = 1; j < cols - 1; ++j)
        {
            temp_grid[i][j] = grid[i][j] + tau * (-((u[i][j] + fabs(u[i][j])) / 2.0) * (grid[i][j] - grid[i - 1][j]) / step - ((u[i][j] - fabs(u[i][j])) / 2.0) * (grid[i + 1][j] - grid[i][j]) / step - ((v[i][j] + fabs(v[i][j])) / 2.0) * (grid[i][j] - grid[i][j - 1]) / step - ((v[i][j] - fabs(v[i][j])) / 2.0) * (grid[i][j + 1] - grid[i][j]) / step + kinematic_viscosity * ((grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1] - 4 * grid[i][j]) / (step * step)));
        }
    }

    copy_vector(grid, temp_grid);
}

void ParallelCalculation::update_velocity(std::vector<std::vector<double>> &grid, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &pressure, double step, double tau, double kinematic_viscosity, double density)
{
    int rows = grid.size();
    int cols = grid[0].size();

    std::vector<std::vector<double>> temp_grid(rows, std::vector<double>(cols, 0.0));

#pragma omp parallel for collapse(2)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j = 1; j < cols - 1; ++j)
        {
            temp_grid[i][j] = grid[i][j] + tau * (-((u[i][j] + std::abs(u[i][j])) / 2 * (grid[i][j] - grid[i - 1][j]) / step + (u[i][j] - std::abs(u[i][j])) / 2 * (grid[i + 1][j] - grid[i][j]) / step) - ((v[i][j] + std::abs(v[i][j])) / 2 * (grid[i][j] - grid[i][j - 1]) / step + (v[i][j] - std::abs(v[i][j])) / 2 * (grid[i][j + 1] - grid[i][j]) / step) - (pressure[i + 1][j + 1] + pressure[i + 1][j - 1] - pressure[i - 1][j + 1] - pressure[i - 1][j - 1]) / (4 * step * density) + kinematic_viscosity * (grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1] - 4 * grid[i][j]) / (step * step));
        }
    }
}

void ParallelCalculation::update_pressure(std::vector<std::vector<double>> &grid, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, double step, double c, double tau)
{
    int rows = grid.size();
    int cols = grid[0].size();

    std::vector<std::vector<double>> temp_grid(rows, std::vector<double>(cols, 0.0));

#pragma omp parallel for collapse(2)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j = 1; j < cols - 1; ++j)
        {
            temp_grid[i][j] = ((u[i + 1][j + 1] + u[i + 1][j - 1]) - (u[i - 1][j - 1] + u[i - 1][j + 1]) + (v[i - 1][j + 1] + v[i + 1][j + 1]) - (v[i - 1][j - 1] + v[i + 1][j - 1])) / (4 * step);
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j = 1; j < cols - 1; ++j)
        {
            grid[i][j] -= c * tau * temp_grid[i][j];
        }
    }
}