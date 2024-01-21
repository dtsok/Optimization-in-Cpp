#ifndef GA_BINARY_HPP_
#define GA_BINARY_HPP_

#include "../common/common_file.hpp"
#include <cstddef>
#include <limits>

class GA_Binary {
	public:
	size_t maxIterations;
	double acc;
	double real_val;

	double (*function)(const size_t, const double *);

	size_t dim = 0;
	size_t bits_per_num = 0;
	size_t total_bits = 0;
	size_t populationSize = 0;
	double l;
	double h;

	int **P = nullptr;
	double *values = nullptr;
	int *best = nullptr;
	size_t best_index = 0;
	double best_value = std::numeric_limits<double>::max();
	int *worst = nullptr;
	size_t worst_index = 0;
	double worst_value = std::numeric_limits<double>::min();

	void initialize(size_t populationSize);
	void evaluate(int **pop, size_t N);
	double decode(int *point, size_t offset = 0);

	double selection_op = 0.2;
	void selection();
	void roulette_wheel_selection();
	void linear_ranking(double *l_bounds);
	void nonlinear_ranking(double *l_bounds);
	void tournament_selection();

	public:
	GA_Binary(size_t dimensions, double (*func)(const size_t, const double *), double l = -2.5, double h = 2.5);
	~GA_Binary();
	void setParameters(size_t iterations, double acc, double real_val);
};

#endif // GA_BINARY_HPP_
