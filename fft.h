#pragma once
#include <complex>
#include <vector>

class FastFourierTransform{
public:
	void hanning(std::vector<std::complex<float> >* data);
	std::vector<std::complex<float> > transform(std::vector<std::complex<float> > data, bool window=true);
	std::vector<std::complex<float> > inverse_transform(std::vector<std::complex<float> > data);
	std::vector<float> float_transform(std::vector<std::complex<float> > data, bool window=true);
	void to_merge_calculation(int stage, int* out);
	void reverse_bit_sort(std::vector<std::complex<float> >* data);

	FastFourierTransform(int length, int exponent);
	~FastFourierTransform();
private:
	int m_Length;
	int m_Exponent;

	int* m_Stages;
	int* m_Total;
	float* m_Hanning;
	
	std::complex<float>* m_Twiddle;
};