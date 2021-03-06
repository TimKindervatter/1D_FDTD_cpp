#pragma once

#include "Device.h"

class BraggGrating : public Device
{
public:
	BraggGrating(floating_point_t max_frequency);

private:
	void construct_grating_arrays(Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& vec, floating_point_t spacer_region_value, floating_point_t grating_layer_1_value, floating_point_t grating_layer_2_value, int num_periods);
};

