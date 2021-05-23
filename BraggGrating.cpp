#include "BraggGrating.h"

#include <cmath>

BraggGrating::BraggGrating(double max_frequency)
{
	double alternating_layer_1_width = 163e-9;
	double alternating_layer_2_width = 122e-9;
	double spacer_region_width = 100e-9;

	double n1 = 1.5; // SiN
	double n2 = 2.0; // SiO2
	double er1 = std::pow(n1, 2); // SiN is non-magnetic (i.e. mur1 = 1.0, so er1 = n1^2)
	double er2 = std::pow(n2, 2); // SiO2 is non-magnetic (i.e. mur2 = 1.0, so er2 = n2^2)
	int num_periods = 15;
	
	construct_grating_arrays(this->m_layer_widths, spacer_region_width, alternating_layer_1_width, alternating_layer_2_width, num_periods);
	construct_grating_arrays(this->m_layer_permittivities, 1.0, er1, er2, num_periods);
	construct_grating_arrays(this->m_layer_permeabilities, 1.0, 1.0, 1.0, num_periods);
	this->m_boundary_refractive_index = 1.0;

	Device::initialize(max_frequency);
}

void BraggGrating::construct_grating_arrays(std::vector<double>& vec, double spacer_region_value, double grating_layer_1_value, double grating_layer_2_value, int num_periods)
{
	int num_grating_layers = 2 * num_periods; // Each period has one SiN layer and one SiO2 layer

	vec.reserve(num_grating_layers + 2); // All grating layers plus one spacer region on each side
	vec.push_back(spacer_region_value);
	for (int i = 0; i < num_periods; ++i)
	{
		// Pattern is Layer 2, Layer 1, Layer2, Layer1, ...
		vec.push_back(grating_layer_2_value);
		vec.push_back(grating_layer_1_value);
	}
	vec.push_back(spacer_region_value);
}
