#include "Utilities.h"
#include "Device.h"

floating_point_t compute_time_step(const Device& device)
{
	floating_point_t nbc = device.get_boundary_refractive_index();
	floating_point_t dz = device.get_grid_resolution();

	floating_point_t time_step = nbc * dz / (2.0f * GlobalConstants::c);

	return time_step;
}

std::string to_lower(std::string s)
{
	std::transform(std::begin(s), std::end(s), std::begin(s), [](auto character) { return std::tolower(character); });

	return s;
}