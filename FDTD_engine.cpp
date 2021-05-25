#include <string>
#include <cmath>
#include <iostream>

#include "ProblemInstances.h"
#include "Utilities.h"
#include "PythonUtilities.h"

using GlobalConstants::c;
using GlobalConstants::pi;

void update_H(std::vector<double>& Hx, std::vector<double>& Ey, std::vector<double>& mHx, double dz, uint32_t Nz)
{
	for (int nz = 0; nz < Nz - 1; ++nz)
	{
		Hx.at(nz) = Hx.at(nz) + mHx.at(nz) * (Ey.at(nz + 1) - Ey.at(nz)) / dz;
	}
}

void update_E(std::vector<double>& Ey, std::vector<double>& Hx, std::vector<double>& mEy, double dz, uint32_t Nz)
{
	for (int nz = 1; nz < Nz; ++nz)
	{
		Ey.at(nz) = Ey.at(nz) + mEy.at(nz) * (Hx.at(nz) - Hx.at(nz - 1)) / dz;
	}
}

void FDTD_engine()
{
	// Define Problem
	BraggGratingProblemInstance problem_instance{};

	double max_frequency = problem_instance.max_frequency;

	uint32_t Nz = problem_instance.device.get_full_grid_size();
	double dz = problem_instance.device.get_grid_resolution();
	std::vector<double> grid = problem_instance.device.get_grid();

	std::vector<double> epsilon_r = problem_instance.device.get_epsilon_r();
	std::vector<double> mu_r = problem_instance.device.get_mu_r();
	std::vector<double> n = problem_instance.device.get_index_of_refraction();

	double dt = compute_time_step(problem_instance.device);

	int source_location = 1;


	// Compute source parameters
	double tau{};
	double t0{};
	if (to_lower(problem_instance.source_type) == "cw")
	{
		tau = 3 / max_frequency;
		t0 = 3 * tau;
	}
	else if (to_lower(problem_instance.source_type) == "pulse")
	{
		tau = 0.5 / max_frequency;
		t0 = 6 * tau;
	}
	
	// Compute number of time steps
	double max_index_of_refraction = max(n, 1.0f);
	double t_prop = max_index_of_refraction * Nz * dz / c;

	double total_runtime = problem_instance.duration_multiplier * (12 * tau + 5 * t_prop);
	int steps = static_cast<int>(std::ceil(total_runtime / dt));

	// Compute source functions for Ey/Hx mode
	std::vector<double> t = arange(0, steps) * dt;
	double A = std::sqrt(epsilon_r.at(source_location) / mu_r.at(source_location));
	double deltat = n.at(source_location) * dz / (2.0f * c) + dt / 2;

	std::vector<double> Eysrc;
	std::vector<double> Hxsrc;

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
		std::vector<double> Ey_pulse = ((t - t0) / tau);
		Eysrc = exp(-pow(Ey_pulse, 2));

		std::vector<double> Hx_pulse = ((t - t0 + deltat) / tau);
		Hxsrc = -A*exp(-pow(Hx_pulse, 2));
	}

	// Initialize update coefficients
	std::vector<double> mEy = (c * dt) / epsilon_r;
	std::vector<double> mHx = (c * dt) / mu_r;

	std::vector<double> Ey(Nz, 0.0);
	std::vector<double> Hx(Nz, 0.0);

	double h2 = 0.0;
	double h1 = 0.0;
	double e2 = 0.0;
	double e1 = 0.0;

	int num_frequencies = problem_instance.num_frequencies; 
	double time_step = compute_time_step(problem_instance.device);

	int m_num_frequencies{ num_frequencies };
	double m_time_step{ time_step };
	std::vector<double> m_frequencies{ linspace(0.0f, problem_instance.max_frequency, num_frequencies) };
	std::vector<std::complex<double>> m_reflected_fourier(num_frequencies, std::complex<double>{0.0, 0.0});
	std::vector<std::complex<double>> m_transmitted_fourier(num_frequencies, std::complex<double>{0.0, 0.0});
	std::vector<std::complex<double>> m_source_fourier(num_frequencies, std::complex<double>{0.0, 0.0});
	std::vector<double> m_reflectance(num_frequencies, 0.0);
	std::vector<double> m_transmittance(num_frequencies, 0.0);
	std::vector<double> m_conservation_of_energy(num_frequencies, 0.0);

	std::complex<double> imaginary_unit{ 0.0, 1.0 };
	auto m_kernel = exp(imaginary_unit * 2.0 * pi * time_step * m_frequencies);

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

		for (int f = 0; f < m_num_frequencies; ++f)
		{
			m_reflected_fourier.at(f) += pow(m_kernel.at(f), T) * Ey.at(0);
			m_transmitted_fourier.at(f) += pow(m_kernel.at(f), T) * Ey.at(Nz - 1);
			m_source_fourier.at(f) += pow(m_kernel.at(f), T) * Eysrc.at(T);
		}

		m_reflectance = squared_magnitude(m_reflected_fourier / m_source_fourier);
		m_transmittance = squared_magnitude(m_transmitted_fourier / m_source_fourier);
		m_conservation_of_energy = m_reflectance + m_transmittance;

		m_reflected_fourier = m_reflected_fourier * m_time_step;
		m_transmitted_fourier = m_transmitted_fourier * m_time_step;
		m_source_fourier = m_source_fourier * m_time_step;
	}

	return;
}