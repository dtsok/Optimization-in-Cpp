#ifndef COMMON_FILE_HPP_
#define COMMON_FILE_HPP_

#include <cmath>
#include <iostream>
#include <cstddef>
#include <limits>
#include <random>
#include <algorithm>
#include <set>
#include <vector>

double l2_norm(const double *X, const size_t i, const size_t j);

double l2_norm_gen(double *x, double *y, const size_t N);

double ELJ(const size_t N, const double *X);

class RandomGenerator {
	public:
	static double generateDouble(double l, double h);
	static int generateInt(double l, double h);

	private:
	static std::default_random_engine generator;
};

// double ELJ(const int N, const double *X, const double epsilon = 1, const double sigma = 1);

#endif // COMMON_FILE_HPP_
