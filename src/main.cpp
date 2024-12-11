#include "ConsoleLogger.h"
#include "ModelFactory.h"
#include "SolverFactory.h"
#include "FluidDynamic.h"
#include <clocale>

int main()
{
    setlocale(LC_ALL, "ru");

    std::shared_ptr<ILogger> logger = std::make_shared<ConsoleLogger>();

    double 
        density = 1000.0,
        kinematic_viscosity = 0.009;

    double 
        Lx = 15.0,
        Ly = 10.0,
        step = 1.0;

    double 
        tmax = 1.5,
        tau = 0.0001,
        u_max = 2.0;

    while (true)
    {
        std::cout << "Methods:\n"
                  << "vorticity-stream\n"
                  << "velocity-pressure\n"
                  << "q\\quit - exit\n";
        std::string method;
        std::cin >> method;
        auto model = ModelFactory::createModel(method);

        std::cout << "Solvers:\n"
                  << "sequential\n"
                  << "parallel\n"
                  << "q\\quit - exit\n";
        std::string calculator;
        std::cin >> calculator;
        auto solver = SolverFactory::createSolver(calculator);

        if (solver == nullptr || model == nullptr)
            return 1;


        std::string pathname_u = method + calculator + "_result_u.txt";
        std::string pathname_v = method + calculator + "_result_v.txt";

        auto solution = FlyidDynamic(std::move(model), std::move(solver), logger);

        logger->log("start init", LogLevel::WARNING);
        solution.set_dim_params(Lx, Ly, step);
        solution.set_model_params(u_max, density, kinematic_viscosity);
        solution.set_time_params(tmax, tau);

        logger->log("Solve", LogLevel::WARNING);
        solution.start();

        logger->log("print resul" + pathname_u + " " + pathname_v, LogLevel::WARNING);
        solution.print_results(pathname_u, pathname_v);
    }

    return 0;
}