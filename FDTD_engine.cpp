#include <string>
#include <cmath>
#include <iostream>

#include "ProblemInstances.h"
#include "Utilities.h"
#include "PythonUtilities.h"

using GlobalConstants::c;
using GlobalConstants::pi;

inline void update_H(std::vector<floating_point_t>& Hx, std::vector<floating_point_t>& Ey, std::vector<floating_point_t>& mHx, floating_point_t dz, uint32_t Nz)
{
	for (size_t nz = 0; nz < Nz - 1; ++nz)
	{
		Hx.at(nz) = Hx.at(nz) + mHx.at(nz) * (Ey.at(nz + 1) - Ey.at(nz)) / dz;
	}
}

inline void update_E(std::vector<floating_point_t>& Ey, std::vector<floating_point_t>& Hx, std::vector<floating_point_t>& mEy, floating_point_t dz, uint32_t Nz)
{
	for (size_t nz = 1; nz < Nz; ++nz)
	{
		Ey.at(nz) = Ey.at(nz) + mEy.at(nz) * (Hx.at(nz) - Hx.at(nz - 1)) / dz;
	}
}

void FDTD_engine()
{
	// Define Problem
	BraggGratingProblemInstance problem_instance{};

	floating_point_t max_frequency = problem_instance.max_frequency;

	uint32_t Nz = problem_instance.device.get_full_grid_size();
	floating_point_t dz = problem_instance.device.get_grid_resolution();
	std::vector<floating_point_t> grid = problem_instance.device.get_grid();

	std::vector<floating_point_t> epsilon_r = problem_instance.device.get_epsilon_r();
	std::vector<floating_point_t> mu_r = problem_instance.device.get_mu_r();
	std::vector<floating_point_t> n = problem_instance.device.get_index_of_refraction();

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
	floating_point_t max_index_of_refraction = max(n, 1.0f);
	floating_point_t t_prop = max_index_of_refraction * Nz * dz / c;

	floating_point_t total_runtime = problem_instance.duration_multiplier * (12 * tau + 5 * t_prop);
	int steps = static_cast<int>(std::ceil(total_runtime / dt));

	// Compute source functions for Ey/Hx mode
	std::vector<floating_point_t> t = arange(0, steps) * dt;
	floating_point_t A = std::sqrt(epsilon_r.at(source_location) / mu_r.at(source_location));
	floating_point_t deltat = n.at(source_location) * dz / (2.0f * c) + dt / 2;

	std::vector<floating_point_t> Eysrc;
	std::vector<floating_point_t> Hxsrc;

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
		std::vector<floating_point_t> Ey_pulse = ((t - t0) / tau);
		Eysrc = exp(-pow(Ey_pulse, 2));

		std::vector<floating_point_t> Hx_pulse = ((t - t0 + deltat) / tau);
		Hxsrc = -A*exp(-pow(Hx_pulse, 2));
	}

	// Initialize update coefficients
	std::vector<floating_point_t> mEy = (c * dt) / epsilon_r;
	std::vector<floating_point_t> mHx = (c * dt) / mu_r;

	std::vector<floating_point_t> Ey(Nz, 0.0);
	std::vector<floating_point_t> Hx(Nz, 0.0);

	floating_point_t h2 = 0.0;
	floating_point_t h1 = 0.0;
	floating_point_t e2 = 0.0;
	floating_point_t e1 = 0.0;

	int num_frequencies = problem_instance.num_frequencies; 
	floating_point_t time_step = compute_time_step(problem_instance.device);

	int m_num_frequencies{ num_frequencies };
	floating_point_t m_time_step{ time_step };
	std::vector<floating_point_t> m_frequencies{ linspace(0.0f, problem_instance.max_frequency, num_frequencies) };
	std::vector<std::complex<floating_point_t>> m_reflected_fourier(num_frequencies, std::complex<floating_point_t>{0.0, 0.0});
	std::vector<std::complex<floating_point_t>> m_transmitted_fourier(num_frequencies, std::complex<floating_point_t>{0.0, 0.0});
	std::vector<std::complex<floating_point_t>> m_source_fourier(num_frequencies, std::complex<floating_point_t>{0.0, 0.0});
	std::vector<floating_point_t> m_reflectance(num_frequencies, 0.0);
	std::vector<floating_point_t> m_transmittance(num_frequencies, 0.0);
	std::vector<floating_point_t> m_conservation_of_energy(num_frequencies, 0.0);

	std::complex<floating_point_t> imaginary_unit{ 0.0, 1.0 };
	auto two_pi_i =  2.0f * pi<floating_point_t> * imaginary_unit;
	auto m_kernel = exp(two_pi_i * time_step * m_frequencies);

	for (int T = 0; T < steps; ++T)
	{
#ifdef _DEBUG
		std::cout << "Starting iteration " << T << " of " << steps << ".\n";
#endif

		// Record H at boundary
		h2 = h1;
		h1 = Hx.at(0);

		// Update H from E
		update_H(Hx, Ey, mHx, dz, Nz);
		Hx.at(Nz - 1) = Hx.at(Nz - 1) + mHx.at(Nz - 1) * (e2 - Ey.at(Nz - 1)) / dz;

		// Handle H-field source
		Hx.at(source_location - 1) = Hx.at(source_location - 1) - (mHx.at(source_location - 1) / dz) * Eysrc.at(T);

		// Record E at boundary
		e2 = e1;
		e1 = Ey.at(Nz - 1);

		// Update E from H
		Ey.at(0) = Ey.at(0) + mEy.at(0) * (Hx.at(0) - h2) / dz;
		update_E(Ey, Hx, mEy, dz, Nz);

		// Handle E-field source
		Ey.at(source_location) = Ey.at(source_location) - (mEy.at(source_location) / dz) * Hxsrc.at(T);


		auto m_kernel_pow = pow(m_kernel, T);
		auto Ey0 = Ey[0];
		auto EyNz = Ey[Nz - 1];
		auto EysrcT = Eysrc[T];

		for (int f = 0; f < m_num_frequencies; ++f)
		{
			m_reflected_fourier[f] += m_kernel_pow[f] * Ey0;
			m_transmitted_fourier[f] += m_kernel_pow[f] * EyNz;
			m_source_fourier[f] += m_kernel_pow[f] * EysrcT;
		}

		auto reflected_fraction = m_reflected_fourier / m_source_fourier;
		m_reflectance = squared_magnitude(reflected_fraction);

		auto transmitted_fraction = m_transmitted_fourier / m_source_fourier;
		m_transmittance = squared_magnitude(transmitted_fraction);
		m_conservation_of_energy = m_reflectance + m_transmittance;

		m_reflected_fourier = m_reflected_fourier * m_time_step;
		m_transmitted_fourier = m_transmitted_fourier * m_time_step;
		m_source_fourier = m_source_fourier * m_time_step;
	}

	return;
}