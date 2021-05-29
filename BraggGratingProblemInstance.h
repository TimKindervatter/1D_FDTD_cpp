#pragma once

#include <string>

#include "BraggGrating.h"
#include "ProblemInstance.h"

struct BraggGratingProblemInstance : public ProblemInstance
{
	floating_point_t max_frequency{ 1.5e15f }; // Hz
	BraggGrating device = BraggGrating(max_frequency);
	int num_frequencies{ 500 };
	int plot_update_interval{ 20 };
	int duration_multiplier{ 3 };
	std::string axis_scaling{ "Logarithmic" };
	floating_point_t ymin{ -2 };
	floating_point_t ymax{ 2 };
	std::string source_type{ "Pulse" };
};