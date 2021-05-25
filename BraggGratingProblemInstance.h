#pragma once

#include <string>

#include "BraggGrating.h"
#include "FourierTransform.h"

struct BraggGratingProblemInstance
{
	double max_frequency{ 1.5e15f }; // Hz
	BraggGrating device = BraggGrating(max_frequency);
	int num_frequencies{ 500 };
	int plot_update_interval{ 20 };
	int duration_multiplier{ 3 };
	std::string axis_scaling{ "Logarithmic" };
	double ymin{ -2 };
	double ymax{ 2 };
	std::string source_type{ "Pulse" };
};