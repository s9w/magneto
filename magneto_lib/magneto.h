#pragma once

#include "export_macro.h"
#include <spdlog/spdlog.h>

namespace magneto {
   CLASS_DECLSPEC std::shared_ptr<spdlog::logger> get_logger(std::vector<spdlog::sink_ptr> sinks = {});
   CLASS_DECLSPEC void start();
}
