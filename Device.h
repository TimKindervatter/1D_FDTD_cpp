#pragma once

#include <vector>
#include <utility>
#include <algorithm>

#include "Utilities.h"

class Device
{
public:
	floating_point_t get_boundary_refractive_index() const { return m_boundary_refractive_index; }
	floating_point_t get_grid_resolution() const { return m_grid_resolution; }
	uint32_t get_full_grid_size() const { return m_full_grid_size; }
	std::vector<floating_point_t> get_grid() const { return m_grid; }
	std::vector<floating_point_t> get_epsilon_r() const { return m_epsilon_r; }
	std::vector<floating_point_t> get_mu_r() const { return m_mu_r; }
	std::vector<floating_point_t> get_index_of_refraction() const { return m_index_of_refraction; }


protected:
	void initialize(floating_point_t max_frequency);

	floating_point_t determine_grid_spacing(floating_point_t max_frequency);
	floating_point_t determine_wavelength_resolution(const std::vector<floating_point_t>& layer_permittivities, const std::vector<floating_point_t>& layer_permeabilities, floating_point_t max_frequency);
	floating_point_t determine_device_resolution(floating_point_t critical_dimension);
	floating_point_t snap_grid(floating_point_t critical_dimension, floating_point_t delta_lambda, floating_point_t delta_d);
	uint32_t round_cells_up(floating_point_t unrounded_number_of_cells);
	std::vector<uint32_t> compute_layer_sizes();
	uint32_t compute_grid_size(const std::vector<uint32_t>& layer_sizes);
	void generate_material_grid(std::vector<floating_point_t>& epsilon_r, std::vector<floating_point_t>& mu_r);
	void generate_grid_1D(std::vector<floating_point_t>& grid);
	std::pair<size_t, size_t> compute_layer_start_and_end_indices(const std::vector<uint32_t> layer_sizes, int i);

	std::vector<floating_point_t> m_layer_widths;
	std::vector<floating_point_t> m_layer_permittivities;
	std::vector<floating_point_t> m_layer_permeabilities;
	std::vector<floating_point_t> m_layer_refractive_indices;
	floating_point_t m_grid_resolution;
	std::vector<uint32_t> m_layer_sizes;
	uint32_t m_full_grid_size;
	std::vector<floating_point_t> m_epsilon_r;
	std::vector<floating_point_t> m_mu_r;
	std::vector<floating_point_t> m_index_of_refraction;
	std::vector<floating_point_t> m_grid;
	floating_point_t m_boundary_refractive_index;
};

