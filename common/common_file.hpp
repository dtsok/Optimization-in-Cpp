#ifndef COMMON_FILE_HPP_
#define COMMON_FILE_HPP_

#include <cmath>
#include <iostream>
#include <random>

float l2_norm(const float *X, const size_t i, const size_t j);

float l2_norm_gen(float *x, float *y, const size_t N);

float ELJ(const size_t N, const float *X);

class RandomGenerator {
	public:
	static float generateFloat(float l, float h);

	private:
	static std::default_random_engine generator;
};

// float ELJ(const int N, const float *X, const float epsilon = 1, const float sigma = 1);

#endif // COMMON_FILE_HPP_