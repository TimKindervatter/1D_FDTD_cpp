#include "Device.h"
#include "Utilities.h"
#include "PythonUtilities.h"

#include <numeric>
#include <cassert>
#include <cmath>

using GlobalConstants::c;

void Device::initialize(floating_point_t max_frequency)
{
	// Base class initialization can only be done after derived class constructor has finished. This function assigns allows deferred "initialization" (members are actually assigned, not initialized)
	// of base class once derived class constructor is finished.
	
	m_layer_refractive_indices = sqrt(m_layer_permittivities * m_layer_permeabilities);
	m_grid_resolution = determine_grid_spacing(max_frequency);
	m_layer_sizes = compute_layer_sizes();
	m_full_grid_size = compute_grid_size(m_layer_sizes);
	generate_material_grid(m_epsilon_r, m_mu_r);
	m_index_of_refraction = sqrt(m_epsilon_r * m_mu_r);
	generate_grid_1D(m_grid);
}

floating_point_t Device::determine_grid_spacing(floating_point_t max_frequency)
{
	floating_point_t delta_lambda = determine_wavelength_resolution(m_layer_permittivities, m_layer_permeabilities, max_frequency);
	auto critical_dimension_iter = std::min_element(std::begin(m_layer_widths), std::end(m_layer_widths));

	floating_point_t critical_dimension{};
	if (critical_dimension_iter != std::end(m_layer_widths))
	{
		critical_dimension = *critical_dimension_iter;
	}

	floating_point_t delta_d = determine_device_resolution(critical_dimension);
	floating_point_t grid_step_size = snap_grid(critical_dimension, delta_lambda, delta_d);

	return grid_step_size;
}

floating_point_t Device::determine_wavelength_resolution(const std::vector<floating_point_t>& layer_permittivities, const std::vector<floating_point_t>& layer_permeabilities, floating_point_t max_frequency)
{
	assert(layer_permittivities.size() == layer_permeabilities.size());

	uint32_t wavelength_resolution = 20;

	std::vector<floating_point_t> layer_indices_of_refraction = sqrt(layer_permittivities * layer_permeabilities);

	floating_point_t max_index_of_refraction = max(layer_indices_of_refraction, 1.0f);

	floating_point_t lambda_min = c / (max_frequency * max_index_of_refraction);
	floating_point_t delta_lambda = lambda_min / wavelength_resolution;

	return delta_lambda;
}

floating_point_t Device::determine_device_resolution(floating_point_t critical_dimension)
{
	floating_point_t device_resolution = 4; // Number of points to resolve device geometry with
	floating_point_t delta_d = critical_dimension / device_resolution;

	return delta_d;
}

floating_point_t Device::snap_grid(floating_point_t critical_dimension, floating_point_t delta_lambda, floating_point_t delta_d)
{
	floating_point_t grid_step_size_unsnapped = std::min(delta_lambda, delta_d);
	floating_point_t grid_size_unsnapped = critical_dimension / grid_step_size_unsnapped;
	uint32_t device_size = round_cells_up(grid_size_unsnapped);
	floating_point_t grid_step_size = critical_dimension / device_size;

	return grid_step_size;
}

uint32_t Device::round_cells_up(floating_point_t unrounded_number_of_cells)
{
	return static_cast<uint32_t>(std::ceil(unrounded_number_of_cells));
}

std::vector<uint32_t> Device::compute_layer_sizes()
{
	std::vector<uint32_t> layer_sizes;
	layer_sizes.reserve(m_layer_widths.size());

	for (auto& layer_width : m_layer_widths)
	{
		layer_sizes.emplace_back(static_cast<uint32_t>(std::ceil(layer_width / m_grid_resolution)));
	}

	return layer_sizes;
}

uint32_t Device::compute_grid_size(const std::vector<uint32_t>& layer_sizes)
{
	uint32_t num_reflection_cells = 1;
	uint32_t num_source_cells = 1;
	uint32_t num_transmission_cells = 1;

	uint32_t num_layer_cells = std::accumulate(std::begin(layer_sizes), std::end(layer_sizes), 0u);

	uint32_t full_grid_size = num_reflection_cells + num_source_cells + num_layer_cells + num_transmission_cells;

	return full_grid_size;
}

void Device::generate_material_grid(std::vector<floating_point_t>& epsilon_r, std::vector<floating_point_t>& mu_r)
{
	assert(m_layer_permittivities.size() == m_layer_sizes.size());
	assert(m_layer_permeabilities.size() == m_layer_sizes.size());

	epsilon_r.resize(m_full_grid_size);
	std::fill(std::begin(epsilon_r), std::end(epsilon_r), 1.0f);

	mu_r.resize(m_full_grid_size);
	std::fill(std::begin(mu_r), std::end(mu_r), 1.0f);

	for (size_t i = 0; i < m_layer_sizes.size(); ++i)
	{
		auto [layer_start_index, layer_end_index] = compute_layer_start_and_end_indices(m_layer_sizes, i);

		std::fill(std::begin(epsilon_r) + layer_start_index, std::begin(epsilon_r) + layer_end_index, m_layer_permittivities.at(i));
		std::fill(std::begin(mu_r) + layer_start_index, std::begin(mu_r) + layer_end_index, m_layer_permeabilities.at(i));
	}
}

void Device::generate_grid_1D(std::vector<floating_point_t>& grid)
{
	floating_point_t zmin = 0.0;
	floating_point_t zmax = m_full_grid_size * m_grid_resolution;

	grid = linspace(zmin, zmax, m_full_grid_size);
}

std::pair<size_t, size_t> Device::compute_layer_start_and_end_indices(const std::vector<uint32_t> layer_sizes, int i)
{
	uint32_t offset = 2; // num_reflection_cells + num_source_cells
	uint32_t layer_start_index = offset + std::accumulate(std::begin(m_layer_sizes), std::begin(m_layer_sizes) + i, 0u);
	uint32_t layer_end_index = offset + std::accumulate(std::begin(m_layer_sizes), std::begin(m_layer_sizes) + i + 1, 0u);

	return std::make_pair(layer_start_index, layer_end_index);
}
