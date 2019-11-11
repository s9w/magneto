#pragma once

#include <spdlog/spdlog.h>
#include "export_macro.h"
#include <memory>


namespace magneto {
   CLASS_DECLSPEC std::shared_ptr<spdlog::logger> get_logger();
}
