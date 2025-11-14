//
// Created by Patrick Flege on 13/11/2025.
//


#include "config/config.hpp"
#include <memory>
#include "cli_caller.hpp"

#include "annotate.hpp"
#include "build.hpp"
#include "transform_annotation.hpp"

std::unique_ptr<mtg::cli::Config> CommandLineInterfaceCaller::get_config(int argc,char** cmd_arguments) {
    return make_unique<mtg::cli::Config>(argc,cmd_arguments);
}

void CommandLineInterfaceCaller::call_flow(const std::unique_ptr<mtg::cli::Config>& config) {
        switch (config->identity) {
            case mtg::cli::Config::BUILD:
                mtg::cli::build_graph(config.get());

            case mtg::cli::Config::ANNOTATE:
                mtg::cli::annotate_graph(config.get());

            case mtg::cli::Config::TRANSFORM_ANNOTATION:
                mtg::cli::transform_annotation(config.get());
            default:
                std::cerr << "Unknown config type" << std::endl;
        }
}
