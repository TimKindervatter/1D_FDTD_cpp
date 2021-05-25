#pragma once

#include <string>

#include "Slab.h"
#include "FourierTransform.h"

struct SlabProblemInstance
{
	double max_frequency{1e9f}; // Hz
	Slab device = Slab(max_frequency);
	int num_frequencies{100};
	int plot_update_interval{20};
	int duration_multiplier{1};
	std::string axis_scaling{"Linear"};
	double ymin{-2};
	double ymax{2};
	std::string source_type{"Pulse"};
};

