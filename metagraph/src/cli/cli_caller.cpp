//
// Created by Patrick Flege on 13/11/2025.
//


#include "config/config.hpp"
#include <memory>
#include "cli_caller.hpp"

#include "annotate.hpp"
#include "build.hpp"
#include "clean.hpp"
#include "transform_annotation.hpp"

std::unique_ptr<mtg::cli::Config> CommandLineInterfaceCaller::get_config(int argc,char** cmd_arguments) {
    return make_unique<mtg::cli::Config>(argc,cmd_arguments);
}

int CommandLineInterfaceCaller::call_flow(const std::unique_ptr<mtg::cli::Config>& config) {
    std::cout << "Identity is: "<< std::endl;
    std::cout << config->identity << std::endl;

    switch (config->identity) {
        case mtg::cli::Config::BUILD:
            return mtg::cli::build_graph(config.get());

        case mtg::cli::Config::ANNOTATE:
            return mtg::cli::annotate_graph(config.get());

        case mtg::cli::Config::TRANSFORM_ANNOTATION:
            return mtg::cli::transform_annotation(config.get());

        case mtg::cli::Config::CLEAN:
            return mtg::cli::clean_graph(config.get());

        default:
            return 1;
        }
}
