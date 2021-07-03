#include "AntiReflectionLayer.h"

AntiReflectionLayer::AntiReflectionLayer(floating_point_t max_frequency)
{
	floating_point_t device_width = 0.3048f;
	floating_point_t spacer_region_width = device_width / 7.1f; // Results in 10 grid cells
	floating_point_t anti_reflection_width = 1.6779e-2;

	this->m_layer_widths.conservativeResize(3);
	this->m_layer_widths << spacer_region_width, anti_reflection_width, device_width, anti_reflection_width, spacer_region_width; // Layer widths in meters

	this->m_layer_permittivities.conservativeResize(3);
	this->m_layer_permittivities << 1.0, 3.46, 12.0, 3.46, 1.0;

	this->m_layer_permeabilities.conservativeResize(3);
	this->m_layer_permeabilities << 1.0, 1.0, 1.0, 1.0, 1.0;

	this->m_boundary_refractive_index = 1.0;

	Device::initialize(max_frequency);
}