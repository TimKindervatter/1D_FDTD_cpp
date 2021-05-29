#pragma once

#include "Utilities.h"
#include "ProblemInstanceIncludes.h"

template <typename InArrayType, typename OutArrayType, typename Callable>
void parallelize_unary_array_operation(InArrayType& input_array, OutArrayType& output_array, Callable f)
{
#pragma omp parallel
	{
		int num_threads = Eigen::nbThreads();
		int thread_ID = omp_get_thread_num();
		int num_elements = input_array.size();
		int elements_per_thread = num_elements / num_threads;

		if (elements_per_thread * num_threads < num_elements)
			elements_per_thread++;

		int start = thread_ID * elements_per_thread;

		int length = elements_per_thread;
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