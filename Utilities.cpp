#include "Utilities.h"

double compute_time_step(const Device& device)
{
	double nbc = device.get_boundary_refractive_index();
	double dz = device.get_grid_resolution();

	double time_step = nbc * dz / (2.0 * GlobalConstants::c);

	return time_step;
}

std::string to_lower(std::string s)
{
	std::transform(std::begin(s), std::end(s), std::begin(s), [](auto character) { return std::tolower(character); });

	return s;
}