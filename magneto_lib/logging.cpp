#include "logging.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>


std::shared_ptr<spdlog::logger> magneto::get_logger() {
   const std::string logger_name = "magneto_logger";
   auto logger = spdlog::get(logger_name);
   if (!logger) {
      std::vector<spdlog::sink_ptr> sinks;
      sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
      sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("log.txt", true));
      logger = std::make_shared<spdlog::logger>(logger_name,
         std::begin(sinks),
         std::end(sinks));
      spdlog::register_logger(logger);
      logger->set_level(spdlog::level::info);
      spdlog::set_pattern("[%T] %^%v%$");
   }
   return logger;
}
