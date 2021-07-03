#pragma once

#include <string>

#include "AntiReflectionLayer.h"

struct AntiReflectionLayerProblemInstance
{
	floating_point_t max_frequency{ 5e9f }; // Hz
	AntiReflectionLayer device = AntiReflectionLayer(max_frequency);
	int num_frequencies{ 500 };
	int plot_update_interval{ 20 };
	int duration_multiplier{ 1 };
	std::string axis_scaling{ "Linear" };
	floating_point_t ymin{ -2 };
	floating_point_t ymax{ 2 };
	std::string source_type{ "Pulse" };
};
