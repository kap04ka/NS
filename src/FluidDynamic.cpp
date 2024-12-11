#include "FluidDynamic.h"
#include "ILogger.h"
#include <iostream>

FlyidDynamic::FlyidDynamic(std::unique_ptr<IEquationModel> model, std::unique_ptr<ICalculation> calculator, std::shared_ptr<ILogger> logger) : model(std::move(model)), calculator(std::move(calculator)), logger(logger)
{
}

void FlyidDynamic::set_dim_params(double width, double height, double step)
{
    int nx = static_cast<int>(width / step) + 1;
    int ny = static_cast<int>(height / step) + 1;
    std::cout << nx << " " <<  ny << std::endl;

    model->set_dims(nx, ny, step);
    logger->log("set_dims passed", LogLevel::WARNING);
}

void FlyidDynamic::set_time_params(double tmax, double tau)
{
    model->set_time(tmax, tau);
    logger->log("set_time_params passed", LogLevel::WARNING);
}

void FlyidDynamic::set_model_params(double u_max, double density, double kinematic_viscosity)
{
    model->set_params(u_max, density, kinematic_viscosity);
    logger->log("set_model_params passed", LogLevel::WARNING);

}

void FlyidDynamic::start()
{   
    logger->log("solve started", LogLevel::WARNING);
    model->solve(std::move(calculator), logger);
    logger->log("solve passed", LogLevel::WARNING);
}

void FlyidDynamic::print_results(const std::string& pathname_u, const std::string& pathname_v) {
    model->print_res(pathname_u, model->get_u());
    model->print_res(pathname_v, model->get_v());
    logger->log("print_results passed", LogLevel::WARNING);
}