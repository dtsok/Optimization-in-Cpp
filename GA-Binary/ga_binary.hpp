#ifndef GA_BINARY_HPP_
#define GA_BINARY_HPP_

#include "../common/common_file.hpp"
#include <cstddef>
#include <limits>
#include <set>

class GA_Binary {
	public:
	GA_Binary(size_t dimensions, double (*func)(const size_t, const double *), double l = -2.5, double h = 2.5);
	~GA_Binary();
	void setParameters(size_t iterations, double acc, double real_val, size_t popSize);
	void minimize(bool tour); // true : use tournament selection - false : use roulette wheel selection

	private:
	size_t dim = 0;									  // dimension
	double (*function)(const size_t, const double *); // objective function
	double l = 0;									  // low for search-space [l, h]
	double h = 0;									  // high for search-space [l, h]

	size_t maxIterations = 0;
	double acc = 0;
	double real_val = 0; // (global) known minimum

	size_t bits_per_num = 0; // bits needed for rep. a number
	size_t total_bits = 0;	 // bits_per_num*dimension

	size_t populationSize = 0; // same for every operator/phase of algorithm

	int **P = nullptr;			// population
	double *p_values = nullptr; // population's evaluated values (function(p[i]))

	int *best = nullptr; // store best vector (optional)
	// size_t best_index = 0;
	double best_value = std::numeric_limits<double>::max();

	int *worst = nullptr; // store worst vector (optional)
	// size_t worst_index = 0;
	double worst_value = std::numeric_limits<double>::min();

	void initialize(); // initialize P and rest arrays

	double decode(int *point,
				  size_t offset = 0); // auxiliary method - convert binary vector to corresponding real in order to evaluate its value
	void evaluate(int **pop, double *val, size_t N); // given Population "pop" evaluate the values and store them in "val"

	int **S = nullptr;			// selected population
	double *s_values = nullptr; // selected's evaluated values

	// double selection_op = 0.2;
	void roulette_wheel_selection(bool linear);
	void linear_ranking(double *l_bounds);
	void nonlinear_ranking(double *l_bounds);
	void tournament_selection(size_t tourSize);

	void crossoverAndUpdate(const size_t p1, const size_t p2, const std::set<size_t> &ind);
	void crossover(double crop);

	void mutation(double mutop);

	void newPopulation();
};

#endif // GA_BINARY_HPP_
