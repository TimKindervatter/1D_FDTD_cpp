#pragma once

#include <cstdint>
#include <string>

#include "Device.h"

namespace GlobalConstants
{
	constexpr uint32_t c = 299792458;  // Speed of light in m / s
	constexpr double pi = 3.14159265358979323846;
}

double compute_time_step(const Device& device);

std::string to_lower(std::string s);