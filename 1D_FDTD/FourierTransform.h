#pragma once

#include <vector>
#include <complex>
#include "Utilities.h"
#include "PythonUtilities.h"

class FourierTransform
{
public:
	FourierTransform(int Nfreq, double max_frequency, double time_step);

	void update_fourier_transforms(int T, const std::vector<double>& Ey, std::vector<double>& Eysrc, uint32_t Nz);
	void finalize_fourier_transforms();

private:
	int m_num_frequencies;
	double m_time_step;
	std::vector<double> m_frequencies;
	std::vector<std::complex<double>> m_kernel;

	std::vector<std::complex<double>> m_reflected_fourier;
	std::vector<std::complex<double>> m_transmitted_fourier;
	std::vector<std::complex<double>> m_source_fourier;

	std::vector<double> m_reflectance;
	std::vector<double> m_transmittance;
	std::vector<double> m_conservation_of_energy;
};

