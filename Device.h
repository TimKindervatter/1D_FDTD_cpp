#pragma once

#include <vector>
#include <utility>
#include <algorithm>

class Device
{
public:
	double get_boundary_refractive_index() const { return m_boundary_refractive_index; }
	double get_grid_resolution() const { return m_grid_resolution; }
	uint32_t get_full_grid_size() const { return m_full_grid_size; }
	std::vector<double> get_grid() const { return m_grid; }
	std::vector<double> get_epsilon_r() const { return m_epsilon_r; }
	std::vector<double> get_mu_r() const { return m_mu_r; }
	std::vector<double> get_index_of_refraction() const { return m_index_of_refraction; }


protected:
	void initialize(double max_frequency);

	double determine_grid_spacing(double max_frequency);
	double determine_wavelength_resolution(const std::vector<double>& layer_permittivities, const std::vector<double>& layer_permeabilities, double max_frequency);
	double determine_device_resolution(double critical_dimension);
	double snap_grid(double critical_dimension, double delta_lambda, double delta_d);
	uint32_t round_cells_up(double unrounded_number_of_cells);
	std::vector<uint32_t> compute_layer_sizes();
	uint32_t compute_grid_size(const std::vector<uint32_t>& layer_sizes);
	void generate_material_grid(std::vector<double>& epsilon_r, std::vector<double>& mu_r);
	void generate_grid_1D(std::vector<double>& grid);
	std::pair<size_t, size_t> compute_layer_start_and_end_indices(const std::vector<uint32_t> layer_sizes, int i);

	std::vector<double> m_layer_widths;
	std::vector<double> m_layer_permittivities;
	std::vector<double> m_layer_permeabilities;
	std::vector<double> m_layer_refractive_indices;
	double m_grid_resolution;
	std::vector<uint32_t> m_layer_sizes;
	uint32_t m_full_grid_size;
	std::vector<double> m_epsilon_r;
	std::vector<double> m_mu_r;
	std::vector<double> m_index_of_refraction;
	std::vector<double> m_grid;
	double m_boundary_refractive_index;
};

