#pragma once

#include <cstdint>
#include <string>

class Device;

using floating_point_t = float;

namespace GlobalConstants
{
	constexpr uint32_t c = 299792458;  // Speed of light in m / s

	template<typename T>
	constexpr T pi = T(3.14159265358979323846);
}

floating_point_t compute_time_step(const Device& device);

std::string to_lower(std::string s);