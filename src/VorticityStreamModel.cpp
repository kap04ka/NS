#include "VorticityStreamModel.h"
#include <fstream>
#include <iomanip>
#include <iostream>

VorticityStreamModel::VorticityStreamModel() {}

void VorticityStreamModel::set_params(double u_max, double density, double kinematic_viscosity)
{
    this->u_max = u_max;
    this->density = density;
    this->kinematic_viscosity = kinematic_viscosity;
}

void VorticityStreamModel::set_time(double tmax, double tau)
{
    this->tmax = tmax;
    this->tau = tau;
}

void VorticityStreamModel::set_dims(int nx, int ny, double step)
{
    this->nx = nx;
    this->ny = ny;
    this->step = step;

    std::cout << "In model: "<< nx << " " <<  ny << std::endl;

    u.resize(this->nx, std::vector<double>(this->ny, 0.0));
    v.resize(this->nx, std::vector<double>(this->ny, 0.0));
    stream_function.resize(this->nx, std::vector<double>(this->ny, 0.0));
    stream_function_new.resize(this->nx, std::vector<double>(this->ny, 0.0));
    vorticity.resize(this->nx, std::vector<double>(this->ny, 0.0));
    vorticity_new.resize(this->nx, std::vector<double>(this->ny, 0.0));
    obstruction.resize(this->nx, std::vector<double>(this->ny, 0.0));
    for (int i = 0; i < this->nx; i++)
    {
        obstruction[i][0] = 1;
        obstruction[i][this->ny - 1] = 1;
    }
}

void VorticityStreamModel::init_boundary_conditions()
{   
    // Граничные значения для скорости (параболический профиль на входе и выходе)
    for (int j = 1; j < ny - 1; j++)
    {
        u[0][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
        u[nx - 1][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
    }

    // Граничные условия для функции тока
    for (int i = 0; i < nx; i++)
    {
        stream_function[i][0] = 0.0; // Нижняя граница
    }

    // Левая и правая границы
    for (int j = 1; j < ny - 1; j++)
    {
        stream_function[0][j] = stream_function[0][j - 1] + u[0][j] * step;
        stream_function[nx - 1][j] = stream_function[nx - 1][j - 1] + u[nx - 1][j] * step;
    }

    // Верхняя граница
    stream_function[0][ny - 1] = stream_function[0][ny - 2];
    for (int i = 1; i < nx; i++)
    {
        stream_function[i][ny - 1] = stream_function[i - 1][ny - 1];
    }

    // Расчет потока
    this->Q = stream_function[0][ny - 1] - stream_function[0][0];

    this->epsilon = 0.05;
}

void VorticityStreamModel::update_boundary_conditions()
{
    for (int i = 0; i < nx; ++i)
    {
        // Верхняя граница
        vorticity[i][ny - 1] = vorticity_new[i][ny - 1] = -(stream_function[i][ny - 1] - stream_function[i][ny - 2]) / (step * step);
        // Нижняя граница
        vorticity[i][0] = vorticity_new[i][0] = -(stream_function[i][1] - stream_function[i][0]) / (step * step);
    }

    // Границы вокруг препятствий
    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            if (obstruction[i][j] == 1)
            {
                if (obstruction[i][j - 1] == 0) // Нижняя граница препятствия
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i][j] - stream_function[i][j - 1]) / (step * step);
                if (obstruction[i][j + 1] == 0) // Верхняя граница препятствия
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i][j + 1] - stream_function[i][j]) / (step * step);
                if (obstruction[i - 1][j] == 0) // Левая граница препятствия
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i][j] - stream_function[i - 1][j]) / (step * step);
                if (obstruction[i + 1][j] == 0) // Правая граница препятствия
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i + 1][j] - stream_function[i][j]) / (step * step);
            }
        }
    }
}

void VorticityStreamModel::solve(std::unique_ptr<ICalculation> calculator, std::shared_ptr<ILogger> logger)
{
    double t{};
    logger->log("solver statred!" + std::to_string(tmax), LogLevel::WARNING);
    init_boundary_conditions();
    logger->log("init gone!", LogLevel::WARNING);
    do
    {
        update_boundary_conditions();
        calculator->helmholtz_solve(vorticity, u, v, step, tau, kinematic_viscosity);
        calculator->jacobi(stream_function, vorticity, epsilon, step);
        calculator->four_point_pattern(u, stream_function, step);
        calculator->reverse_four_point_pattern(v, stream_function, step);
        t += tau;
    } while (t <= tau);
    
}

void VorticityStreamModel::print_res(const std::string& pathname, std::vector<std::vector<double>> grid) {
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