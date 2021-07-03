#pragma once

#include <vector>
#include <utility>
#include <algorithm>

#include "Utilities.h"

#include "Eigen\Dense"

class Device
{
public:
	floating_point_t get_boundary_refractive_index() const { return m_boundary_refractive_index; }
	floating_point_t get_grid_resolution() const { return m_grid_resolution; }
	uint32_t get_full_grid_size() const { return m_full_grid_size; }
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> get_grid() const { return m_grid; }
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> get_epsilon_r() const { return m_epsilon_r; }
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> get_mu_r() const { return m_mu_r; }
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> get_index_of_refraction() const { return m_index_of_refraction; }


protected:
	void initialize(floating_point_t max_frequency);

	floating_point_t determine_grid_spacing(floating_point_t max_frequency);
	floating_point_t determine_wavelength_resolution(const Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& layer_permittivities, const Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& layer_permeabilities, floating_point_t max_frequency);
	floating_point_t determine_device_resolution(floating_point_t critical_dimension);
	floating_point_t snap_grid(floating_point_t critical_dimension, floating_point_t delta_lambda, floating_point_t delta_d);
	uint32_t round_cells_up(floating_point_t unrounded_number_of_cells);
	Eigen::Array<uint32_t, 1, Eigen::Dynamic> compute_layer_sizes();
	uint32_t compute_grid_size(const Eigen::Array<uint32_t, 1, Eigen::Dynamic>& layer_sizes);
	void generate_material_grid(Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& epsilon_r, Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& mu_r);
	void generate_grid_1D(Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& grid);
	std::pair<size_t, size_t> compute_layer_start_and_end_indices(const Eigen::Array<uint32_t, 1, Eigen::Dynamic> layer_sizes, int i);

	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_layer_widths;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_layer_permittivities;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_layer_permeabilities;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_layer_refractive_indices;
	floating_point_t m_grid_resolution{};
	Eigen::Array<uint32_t, 1, Eigen::Dynamic> m_layer_sizes;
	uint32_t m_full_grid_size{};
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_epsilon_r;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_mu_r;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_index_of_refraction;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_grid;
	floating_point_t m_boundary_refractive_index;
};

