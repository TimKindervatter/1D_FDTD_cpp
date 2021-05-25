#include "FourierTransform.h"

using GlobalConstants::pi;

FourierTransform::FourierTransform(int Nfreq, double max_frequency, double time_step) :
	m_num_frequencies{ Nfreq },
	m_time_step{ time_step },
	m_frequencies{ linspace(0.0, max_frequency, Nfreq) },
	m_reflected_fourier(Nfreq, std::complex<double>{0.0, 0.0}),
	m_transmitted_fourier(Nfreq, std::complex<double>{0.0, 0.0}),
	m_source_fourier(Nfreq, std::complex<double>{0.0, 0.0}),
	m_reflectance(Nfreq, 0.0),
	m_transmittance(Nfreq, 0.0),
	m_conservation_of_energy(Nfreq, 0.0)
{
	std::complex<double> imaginary_unit{ 0.0, 1.0 };
	m_kernel = exp(imaginary_unit * 2.0 * pi * time_step * m_frequencies);
}

void FourierTransform::update_fourier_transforms(int T, const std::vector<double>& Ey, std::vector<double>& Eysrc, uint32_t Nz)
{
	for (int f = 0; f < m_num_frequencies; ++f)
	{
		m_reflected_fourier.at(f) += pow(m_kernel.at(f), T) * Ey.at(0);
		m_transmitted_fourier.at(f) += pow(m_kernel.at(f), T) * Ey.at(Nz - 1);
		m_source_fourier.at(f) += pow(m_kernel.at(f), T) * Eysrc.at(T);
	}

	m_reflectance = squared_magnitude(m_reflected_fourier / m_source_fourier);
	m_transmittance = squared_magnitude(m_transmitted_fourier / m_source_fourier);
	m_conservation_of_energy = m_reflectance + m_transmittance;
}

void FourierTransform::finalize_fourier_transforms()
{
	m_reflected_fourier = m_reflected_fourier * m_time_step;
	m_transmitted_fourier = m_transmitted_fourier * m_time_step;
	m_source_fourier = m_source_fourier * m_time_step;
}
