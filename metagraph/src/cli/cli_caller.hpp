//
// Created by Patrick Flege on 13/11/2025.
//

#pragma once
#include "cli_caller_interface.hpp"
#include "config/config.hpp"


class CommandLineInterfaceCaller : public AbstractCommandLineInterface{
    public:
        std::unique_ptr<mtg::cli::Config> get_config(int argc, char** cmd_arguments) override;
        void call_flow(const std::unique_ptr<mtg::cli::Config>& config) override;
};
