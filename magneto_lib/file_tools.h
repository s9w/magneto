#pragma once

#include <filesystem>
#include <optional>

namespace magneto {
	std::optional<std::string> get_file_contents(const std::filesystem::path& path);
	void write_string_to_file(const std::filesystem::path& path, const std::string& content);
}
