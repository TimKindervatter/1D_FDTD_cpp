#include "Slab.h"

Slab::Slab(floating_point_t max_frequency)
{
	floating_point_t device_width = 0.3048f;
	floating_point_t spacer_region_width = device_width / 7.1f;

	this->m_layer_widths = { spacer_region_width, device_width, spacer_region_width };
	this->m_layer_permittivities = { 1.0, 6.0, 1.0 };
	this->m_layer_permeabilities = { 1.0, 2.0, 1.0 };
	this->m_boundary_refractive_index = 1.0;

	Device::initialize(max_frequency);
}