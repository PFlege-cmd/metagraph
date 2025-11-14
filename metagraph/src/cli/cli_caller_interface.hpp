//
// Created by Patrick Flege on 13/11/2025.
//

#ifndef CLI_CALLER_INTERFACE_H
#define CLI_CALLER_INTERFACE_H
#include "config/config.hpp"

#endif //CLI_CALLER_INTERFACE_H
#pragma once

class AbstractCommandLineInterface {
    public:
    AbstractCommandLineInterface()= default;
    ~AbstractCommandLineInterface()= default;
    virtual std::unique_ptr<mtg::cli::Config>get_config(int argc,char** cmd_arguments) = 0;
    virtual void call_flow(const unique_ptr<mtg::cli::Config>& config) = 0;
};