#define _USE_MATH_DEFINES
#include "fft.h"
#include <cmath>
#include <iostream>

FastFourierTransform::FastFourierTransform(int length, int exponent){
	m_Length = length;
	m_Exponent = exponent;
	m_Stages = new int[exponent * length / 2];
	m_Twiddle = new std::complex<float>[exponent * length / 2];
	m_Total = new int[length];
	m_Hanning = new float[length];
	int* odd = new int[length / 2];
	int* even = new int[length / 2];
	int* oddcopy = new int[length / 2];
	int* evencopy = new int[length / 2];
	for (int i = 0; i < length; i++){
		m_Total[i] = i;
	}
	for (int x = 0; x < length; x++){
		if (x % 2 == 0){
			even[x / 2] = m_Total[x];
			evencopy[x / 2] = m_Total[x];
		}else{
			odd[x / 2] = m_Total[x];
			oddcopy[x / 2] = m_Total[x];
		}	
	}
	for (int i = 0; i < exponent - 1; i++){
		if (i == 0){
			for (int x = 0; x < length / 2; x++){
				if (x % 2 == 0){
					even[x / 2] = evencopy[x];
					odd[x / 2] = oddcopy[x];
				}else{
					even[x / 2 + length / 4] = evencopy[x];
					odd[x / 2 + length / 4] = oddcopy[x];
				}	
			}
		}else{
			for (int j = 0; j < pow(2.0f, i); j++){
				int section = length / (2 * pow(2.0f, i));
				int starting_point = j * section;
				int ending_point = starting_point + section;
				for (int x = starting_point; x < ending_point; ++x){
					if (x % 2 == 0){
						even[starting_point / 2 + x / 2] = evencopy[x];
						odd[starting_point / 2 + x / 2] = oddcopy[x];
					}else{
						even[starting_point / 2 + x / 2 + section / 2] = evencopy[x];
						odd[starting_point / 2 + x / 2 + section / 2] = oddcopy[x];
					}
				}
			}
		}
		for (int j = 0; j < length / 2; j++){
			oddcopy[j] = odd[j];
			evencopy[j] = even[j];
		}
	}
	
	for (int i = 0; i < length / 2; i++){
		m_Total[i] = even[i];
		m_Total[i + length / 2] = odd[i];
	}

	delete[] odd;
	delete[] even;
	delete[] oddcopy;
	delete[] evencopy;

	int* indices = new int[length / 2];
	std::complex<float> base = std::complex<float>(0.0f, -2.0f);
	for (int stage = 0; stage < exponent; stage++){
		int max_i = pow(2.0f, stage);
		float group_size = max_i * 2;
		to_merge_calculation(stage, indices);
		for (int i = 0; i < length / 2; i++){
			m_Stages[i + stage * length / 2] = m_Total[indices[i]];
			m_Twiddle[i + stage * length / 2] = std::exp(base * (float)M_PI * (float)(i % max_i) / group_size);
		}
	}
	for (int i = 0; i < length; i++){
		m_Hanning[i] = 0.5 * (1 - cos((2 * M_PI * i) / (length  - 1)));
	}
	delete[] indices;
}

FastFourierTransform::~FastFourierTransform(){
	delete[] m_Stages;
	delete[] m_Total;
	delete[] m_Hanning;
	delete[] m_Twiddle;
}

std::vector<std::complex<float> > FastFourierTransform::transform(std::vector<std::complex<float> > data, bool window){
	if (window)
		hanning(&data);
	std::complex<float>* T = new std::complex<float>[m_Length / 2];
	for (int stage = 0; stage < m_Exponent; stage++){
		int index_offset = pow(2.0f, m_Exponent - stage - 1);
		for (int v = 0; v < m_Length / 2; v++){
			T[v] = m_Twiddle[stage * m_Length / 2 + v] * data[m_Stages[stage * m_Length / 2 + v] + index_offset];
		}
		for (int i = 0; i < m_Length / 2; i++){
			int v = m_Stages[stage * m_Length / 2 + i];
			std::complex<float> carry = data[v];
			data[v] = carry + T[i];
 			data[v + index_offset] = carry - T[i];
		}
	}
	reverse_bit_sort(&data);
	delete[] T;
	return data;
}

std::vector<float> FastFourierTransform::float_transform(std::vector<std::complex<float> > data, bool window){
	data = transform(data, window);
	std::vector<float> newdata;
	for (int i = 0; i < m_Length / 2; i++){
		newdata.push_back(std::abs(data[i]));
	}
	return newdata;
}

std::vector<std::complex<float> > FastFourierTransform::inverse_transform(std::vector<std::complex<float > > data){
	for (int i = 0; i < m_Length; ++i){
		data[i] = std::conj(data[i]);
	}
	data = transform(data, false);
	for (int i = 0; i < m_Length; ++i){
		data[i] = std::conj(data[i]);
		data[i] /= m_Length;
	}
	return data;
}

void FastFourierTransform::reverse_bit_sort(std::vector<std::complex<float> >* data){
	std::complex<float>* new_array = new std::complex<float>[m_Length];
	for (int i = 0; i < m_Length; i++){
		new_array[i] = (*data)[m_Total[i]];
	}
	for (int i = 0; i < m_Length; i++){
		(*data)[i] = new_array[i];
	}
	delete[] new_array;
}

void FastFourierTransform::hanning(std::vector<std::complex<float> >* data){
	for (int i = 0; i < m_Length; i++){
		(*data)[i] *= m_Hanning[i];
	}
}

void FastFourierTransform::to_merge_calculation(int stage, int* out){
	int i = 0;
	int j = 0;
	int offset = pow(2.0f, stage);
	while (i < m_Length){
		out[j] = i;
		j++;
		if (i % offset == offset - 1){
			i += offset + 1;
		}else{
			i += 1;
		}
	}
}