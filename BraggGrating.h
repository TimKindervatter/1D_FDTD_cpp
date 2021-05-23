#pragma once

#include "Device.h"

class BraggGrating : public Device
{
public:
	BraggGrating(double max_frequency);

private:
	void construct_grating_arrays(std::vector<double>& vec, double spacer_region_value, double grating_layer_1_value, double grating_layer_2_value, int num_periods);
};

