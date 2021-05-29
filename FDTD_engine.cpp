#include <string>
#include <cmath>
#include <iostream>
#include <thread>

#include "Eigen\Dense"
#include "matplotlibcpp.h"

#include "ProblemInstanceIncludes.h"
#include "Utilities.h"

#include "FDTD_engine.h"

using GlobalConstants::c;
using GlobalConstants::two_pi;

namespace plt = matplotlibcpp;

inline void update_H(Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& Hx, Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& Ey, Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& mHx, floating_point_t dz, uint32_t Nz)
{
	for (size_t nz = 0; nz < Nz - 1; ++nz)
	{
		Hx[nz] = Hx[nz] + mHx[nz] * (Ey[nz + 1] - Ey[nz]) / dz;
	}
}

inline void update_E(Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& Ey, Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& Hx, Eigen::Array<floating_point_t, 1, Eigen::Dynamic>& mEy, floating_point_t dz, uint32_t Nz)
{
	for (size_t nz = 1; nz < Nz; ++nz)
	{
		Ey[nz] = Ey[nz] + mEy[nz] * (Hx[nz] - Hx[nz - 1]) / dz;
	}
}

void plot_fields(int T, Eigen::Array<floating_point_t, 1, -1>& grid, Eigen::Array<floating_point_t, 1, -1>& Ey, Eigen::Array<floating_point_t, 1, -1>& Hx, SlabProblemInstance& problem_instance)
{
	if (T % problem_instance.plot_update_interval == 0)
	{
		plt::clf();
		plt::subplot(2, 1, 1);
		plt::plot(grid, Ey);
		plt::subplot(2, 1, 2);
		plt::plot(grid, Hx);
		plt::xlim(grid[0], grid[grid.size() - 1]);
		plt::ylim(problem_instance.ymin, problem_instance.ymax);
		plt::pause(0.000000000001);
	}
}

void FDTD_engine()
{
	plt::backend("WXAgg");

	// Define Problem
	SlabProblemInstance problem_instance{};

	floating_point_t max_frequency = problem_instance.max_frequency;

	uint32_t Nz = problem_instance.device.get_full_grid_size();
	floating_point_t dz = problem_instance.device.get_grid_resolution();
	auto grid = problem_instance.device.get_grid();

	auto epsilon_r = problem_instance.device.get_epsilon_r();
	auto mu_r = problem_instance.device.get_mu_r();
	auto n = problem_instance.device.get_index_of_refraction();

	floating_point_t dt = compute_time_step(problem_instance.device);

	int source_location = 1;


	// Compute source parameters
	floating_point_t tau{};
	floating_point_t t0{};
	if (to_lower(problem_instance.source_type) == "cw")
	{
		tau = 3 / max_frequency;
		t0 = 3 * tau;
	}
	else if (to_lower(problem_instance.source_type) == "pulse")
	{
		tau = 0.5f / max_frequency;
		t0 = 6 * tau;
	}
	
	// Compute number of time steps
	floating_point_t max_index_of_refraction = n.maxCoeff();
	floating_point_t t_prop = max_index_of_refraction * Nz * dz / c;

	floating_point_t total_runtime = problem_instance.duration_multiplier * (12 * tau + 5 * t_prop);
	int steps = static_cast<int>(std::ceil(total_runtime / dt));

	// Compute source functions for Ey/Hx mode
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> t(steps);
	t = t.setLinSpaced(0, steps) * dt;

	floating_point_t A = std::sqrt(epsilon_r[source_location] / mu_r[source_location]);
	floating_point_t deltat = n[source_location] * dz / (2.0f * c) + dt / 2;

	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> Eysrc;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> Hxsrc;

	if (to_lower(problem_instance.source_type) == "cw")
	{
		//f_op = 1.5e9

		//g_E = np.sin(2 * np.pi * f_op * t) + np.finfo(float).eps
		//g_H = -np.sin(2 * np.pi * f_op * t) + np.finfo(float).eps

		//Eysrc = g_E * np.ones(len(t))
		//Hxsrc = g_H * np.ones(len(t))

		//Eysrc[t < t0] *= np.exp(-((t[t < t0] - t0) / tau) * *2)
		//Hxsrc[t < t0 + deltat] *= A * np.exp(-((t[t < t0 + deltat] - t0 + deltat) / tau) * *2)
	}
	else if (to_lower(problem_instance.source_type) == "pulse")
	{
		Eigen::Array<floating_point_t, 1, Eigen::Dynamic> Ey_pulse = ((t - t0) / tau);
		Eysrc = Eigen::exp(-(Ey_pulse * Ey_pulse));

		Eigen::Array<floating_point_t, 1, Eigen::Dynamic> Hx_pulse = ((t - t0 + deltat) / tau);
		Hxsrc = -A*Eigen::exp(-(Hx_pulse * Hx_pulse));
	}

	// Initialize update coefficients
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> mEy = (c * dt) / epsilon_r;
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> mHx = (c * dt) / mu_r;

	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> Ey(Nz);
	Ey.fill(0.0);
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> Hx(Nz);
	Hx.fill(0.0);

	floating_point_t h2 = 0.0;
	floating_point_t h1 = 0.0;
	floating_point_t e2 = 0.0;
	floating_point_t e1 = 0.0;

	int num_frequencies = problem_instance.num_frequencies; 

	int m_num_frequencies{ num_frequencies };
	floating_point_t time_step{ compute_time_step(problem_instance.device) };
	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_frequencies;
	m_frequencies.conservativeResize(num_frequencies);
	m_frequencies.setLinSpaced(num_frequencies, 0.0f, problem_instance.max_frequency);

	Eigen::Array<std::complex<floating_point_t>, 1, Eigen::Dynamic> m_reflected_fourier(num_frequencies);
	m_reflected_fourier.fill(0.0);

	Eigen::Array<std::complex<floating_point_t>, 1, Eigen::Dynamic> m_transmitted_fourier(num_frequencies);
	m_transmitted_fourier.fill(0.0);

	Eigen::Array<std::complex<floating_point_t>, 1, Eigen::Dynamic> m_source_fourier(num_frequencies);
	m_source_fourier.fill(0.0);

	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_reflectance(num_frequencies);
	m_reflectance.fill(0.0);

	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_transmittance(num_frequencies);
	m_transmittance.fill(0.0);

	Eigen::Array<floating_point_t, 1, Eigen::Dynamic> m_conservation_of_energy(num_frequencies);
	m_conservation_of_energy.fill(0.0);

	std::complex<floating_point_t> imaginary_unit{ 0.0, 1.0 };
	Eigen::Array<std::complex<floating_point_t>, 1, Eigen::Dynamic> two_pi_i_f = imaginary_unit * two_pi<floating_point_t> * time_step * m_frequencies;

	for (int T = 0; T < steps; ++T)
	{
#ifdef _DEBUG
		std::cout << "Starting iteration " << T << " of " << steps << ".\n";
#endif

		// Record H at boundary
		h2 = h1;
		h1 = Hx[0];

		// Update H from E
		update_H(Hx, Ey, mHx, dz, Nz);
		Hx[Nz - 1] = Hx[Nz - 1] + mHx[Nz - 1] * (e2 - Ey[Nz - 1]) / dz;

		// Handle H-field source
		Hx[source_location - 1] = Hx[source_location - 1] - (mHx[source_location - 1] / dz) * Eysrc[T];

		// Record E at boundary
		e2 = e1;
		e1 = Ey[Nz - 1];

		// Update E from H
		Ey[0] = Ey[0] + mEy[0] * (Hx[0] - h2) / dz;
		update_E(Ey, Hx, mEy, dz, Nz);

		// Handle E-field source
		Ey[source_location] = Ey[source_location] - (mEy[source_location] / dz) * Hxsrc[T];

		// Stupidly, this is faster than .pow(T) because Eigen vectorizes *, and .exp(), but not .pow()
		Eigen::Array<std::complex<floating_point_t>, 1, Eigen::Dynamic> two_pi_i_f_T = two_pi_i_f * T;
		Eigen::Array<std::complex<floating_point_t>, 1, Eigen::Dynamic> m_kernel_pow(two_pi_i_f_T.size());

		parallelize_unary_array_operation(two_pi_i_f_T, m_kernel_pow, [](auto array) { return Eigen::exp(array); });

		auto Ey0 = Ey[0];
		auto EyNz = Ey[Nz - 1];
		auto EysrcT = Eysrc[T];

		for (int f = 0; f < m_num_frequencies; ++f)
		{
			m_reflected_fourier[f] += m_kernel_pow[f] * Ey0;
			m_transmitted_fourier[f] += m_kernel_pow[f] * EyNz;
			m_source_fourier[f] += m_kernel_pow[f] * EysrcT;
		}

		Eigen::Matrix<std::complex<floating_point_t>, 1, Eigen::Dynamic> reflected_fraction = m_reflected_fourier / m_source_fourier;
		m_reflectance = reflected_fraction.squaredNorm();

		Eigen::Matrix<std::complex<floating_point_t>, 1, Eigen::Dynamic> transmitted_fraction = m_transmitted_fourier / m_source_fourier;
		m_transmittance = transmitted_fraction.squaredNorm();
		m_conservation_of_energy = m_reflectance + m_transmittance;

		plot_fields(T, grid, Ey, Hx, problem_instance);
	}

	m_reflected_fourier = m_reflected_fourier * time_step;
	m_transmitted_fourier = m_transmitted_fourier * time_step;
	m_source_fourier = m_source_fourier * time_step;

#ifdef DEBUG
	std::cout << m_reflected_fourier << "\n\n";
	std::cout << m_transmitted_fourier << "\n\n";
	std::cout << m_source_fourier << "\n\n";
#endif
	return;
}