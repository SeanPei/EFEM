#include "CLI/CLI11.hpp"
#include "simulation/Simulator.h"

int main(int argc, char** argv)
{
    CLI::App app { "Efficient FEM" };
    std::string config_file;
    app.add_option("-f", config_file, "Config file");

    CLI11_PARSE(app, argc, argv);
    spdlog::info("config_file={}", config_file);

    using namespace EFEM;
    auto simulator = std::make_shared<Simulator>();
    simulator->initial(config_file);
    simulator->start();
    return 0;
}