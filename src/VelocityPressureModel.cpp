#include "VelocityPressureModel.h"
#include <fstream>
#include <iomanip>
#include <iostream>
VelocityPressureModel::VelocityPressureModel() {}

void VelocityPressureModel::set_params(double u_max, double density, double kinematic_viscosity)
{
    this->u_max = u_max;
    this->density = density;
    this->kinematic_viscosity = kinematic_viscosity;
    this->c = 1;
}

void VelocityPressureModel::set_time(double tmax, double tau)
{
    this->tau = tau;
    this->tmax = tmax;
}

void VelocityPressureModel::set_dims(int nx, int ny, double step)
{
    this->nx = nx;
    this->ny = ny;
    this->step = step;

    u.resize(this->nx, std::vector<double>(this->ny, 0.0));
    v.resize(this->nx, std::vector<double>(this->ny, 0.0));
    u_new.resize(this->nx, std::vector<double>(this->ny, 0.0));
    v_new.resize(this->nx, std::vector<double>(this->ny, 0.0));
    pressure.resize(this->nx, std::vector<double>(this->ny, 0.0));
    div_velocity.resize(this->nx, std::vector<double>(this->ny, 0.0));
}

void VelocityPressureModel::init_boundary_conditions()
{
    for (int j = 1; j < ny - 1; j++)
    {
        u[0][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
        u[nx - 1][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
    }
}

void VelocityPressureModel::update_boundary_conditions()
{
    // Грани
    for (int i = 1; i < nx - 1; ++i)
    {
        pressure[i][0] = pressure[i][1];
        pressure[i][ny - 1] = pressure[i][ny - 2];
    }
    for (int j = 1; j < ny - 1; ++j)
    {
        pressure[0][j] = 2 * pressure[1][j] - pressure[2][j];
        pressure[nx - 1][j] = 2 * pressure[nx - 2][j] - pressure[nx - 3][j];
    }
}

void VelocityPressureModel::solve(std::unique_ptr<ICalculation> calculator, std::shared_ptr<ILogger> logger)
{
    double t{};
    logger->log("solver statred!" + std::to_string(tmax), LogLevel::WARNING);
    init_boundary_conditions();
    logger->log("init gone!", LogLevel::WARNING);
    do
    {
        calculator->update_pressure(pressure, u, v, step, c, tau);
        calculator->update_velocity(u, u, v, pressure, step, tau, kinematic_viscosity, density);
        calculator->update_velocity(v, u, v, pressure, step, tau, kinematic_viscosity, density);
        t += tau;
    } while (t <= tmax);
}

void VelocityPressureModel::print_res(const std::string& pathname, std::vector<std::vector<double>> grid) {
    std::ofstream out;
    out.open(pathname);
    for(int j = ny - 1; j >= 0; --j) {
        for (int i = 0; i < nx; ++i) {
            out << std::fixed << std::setprecision(6) << grid[i][j] << "\t\t";
        }
        out << std::endl;
    }
    out.close();
}