#pragma once

#include "Utilities.h"
#include "ProblemInstanceIncludes.h"

//#define PLOT

// See here for info on this openMP parallelization code: https://stackoverflow.com/questions/41845724/eigen-openmp-no-parallelisation-due-to-false-sharing-and-thread-overhead

template <typename InArrayType, typename OutArrayType, typename Callable>
void parallelize_unary_array_operation(InArrayType& input_array, OutArrayType& output_array, Callable f)
{
#pragma omp parallel
	{
		int num_threads = Eigen::nbThreads();
		int thread_ID = omp_get_thread_num();
		auto num_elements = input_array.size();
		auto elements_per_thread = num_elements / num_threads;

		if (elements_per_thread * num_threads < num_elements)
			elements_per_thread++;

		auto start = thread_ID * elements_per_thread;

		auto length = elements_per_thread;
		if (thread_ID + 1 == num_threads)
			length = num_elements - start;

		if (start < num_elements)
		{
			output_array.segment(start, length) = f(input_array.segment(start, length));
		}
	}
}

template <typename InArray1Type, typename InArray2Type, typename OutArrayType, typename Callable>
void parallelize_binary_array_operation(InArray1Type& input_array1, InArray2Type& input_array2, OutArrayType& output_array, Callable f)
{
#pragma omp parallel
	{
		int num_threads = Eigen::nbThreads();
		int thread_ID = omp_get_thread_num();
		int num_elements = input_array1.size();
		int elements_per_thread = num_elements / num_threads;

		if (elements_per_thread * num_threads < num_elements)
			elements_per_thread++;

		int start = thread_ID * elements_per_thread;

		int length = elements_per_thread;
		if (thread_ID + 1 == num_threads)
			length = num_elements - start;

		if (start < num_elements)
		{
			output_array.segment(start, length) = f(input_array1.segment(start, length), input_array2.segment(start, length));
		}
	}
}

void FDTD_engine();

template <typename ProblemInstance>
void plot_fields(int T,
	Eigen::Array<floating_point_t, 1, -1>& grid,
	Eigen::Array<floating_point_t, 1, -1>& Ey,
	Eigen::Array<floating_point_t, 1, -1>& Hx,
	ProblemInstance& problem_instance,
	Eigen::Array<floating_point_t, 1, -1>& frequencies,
	Eigen::Array<floating_point_t, 1, -1>& reflectance,
	Eigen::Array<floating_point_t, 1, -1>& transmittance,
	Eigen::Array<floating_point_t, 1, -1>& conservation_of_energy)
{
	if (T % problem_instance.plot_update_interval == 0)
	{


		plt::clf();

		plt::subplot(2, 1, 1);
		plt::plot(grid, Ey);
		plt::plot(grid, Hx);
		plt::xlim(grid[0], grid[grid.size() - 1]);
		plt::ylim(problem_instance.ymin, problem_instance.ymax);

		if (to_lower(problem_instance.axis_scaling) == "linear")
		{
			plt::subplot(2, 1, 2);
			plt::plot(frequencies, reflectance, std::map<std::string, std::string>{ {"label", "Reflectance"}});
			plt::plot(frequencies, transmittance, std::map<std::string, std::string>{ {"label", "Transmittance"}});
			plt::plot(frequencies, conservation_of_energy, std::map<std::string, std::string>{ {"label", "Conservation"}});
			plt::xlim(static_cast<floating_point_t>(0.0), problem_instance.max_frequency);
			plt::ylim(0.0, 1.5);
			plt::legend();
		}
		else if (to_lower(problem_instance.axis_scaling) == "logarithmic")
		{
			Eigen::Array<floating_point_t, 1, -1> log_reflectance = 10 * reflectance.log10();
			Eigen::Array<floating_point_t, 1, -1> log_transmittance = 10 * transmittance.log10();
			Eigen::Array<floating_point_t, 1, -1> log_conservation_of_energy = 10 * conservation_of_energy;

			plt::subplot(2, 1, 2);
			plt::plot(frequencies, log_reflectance, std::map<std::string, std::string>{ {"label", "Reflectance"}});
			plt::plot(frequencies, log_transmittance, std::map<std::string, std::string>{ {"label", "Transmittance"}});
			plt::plot(frequencies, log_conservation_of_energy, std::map<std::string, std::string>{ {"label", "Conservation"}});
			plt::legend();
		}

		plt::pause(0.000000000001);
	}
}
