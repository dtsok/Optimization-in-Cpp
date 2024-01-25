#ifndef GA_REAL_HPP_
#define GA_REAL_HPP_

#include "../common/common_file.hpp"
#include <cstddef>
#include <set>

class GA_Real {
	public:
	GA_Real(size_t dimensions, double (*func)(const size_t, const double *), double l = -2.5, double h = 2.5);
	~GA_Real();
	void setParameters(size_t iterations, double acc, double real_val, size_t popSize);
	void minimize(bool tour); // true : use tournament selection - false : use roulette wheel selection

	public:
	size_t dim = 0;									  // dimension
	double (*function)(const size_t, const double *); // objective function
	double l = 0;									  // low for search-space [l, h]
	double h = 0;									  // high for search-space [l, h]

	size_t maxIterations = 0;
	double acc = 0;
	double real_val = 0; // (global) known minimum

	size_t populationSize = 0; // same for every operator/phase of algorithm

	double **P = nullptr;			// population
	double *p_values = nullptr; // population's evaluated values (function(p[i]))

	double *best = nullptr; // store best vector (optional)
	double best_value = std::numeric_limits<double>::max();

	double *worst = nullptr; // store worst vector (optional)
	double worst_value = std::numeric_limits<double>::min();

	void initialize(); // initialize P and rest arrays

	void evaluate(double **pop, double *val, size_t N); // given Population "pop" evaluate the values and store them in "val"

	double **S = nullptr;			// selected population
	double *s_values = nullptr; // selected's evaluated values

	void roulette_wheel_selection(bool linear);
	void linear_ranking(double *l_bounds);
	void nonlinear_ranking(double *l_bounds);
	void tournament_selection(size_t tourSize);

	void crossoverAndUpdate(const size_t p1, const size_t p2, const double delta=0.25);
	void crossover(double crop);

	void mutation(double mutop);

	void newPopulation();
};

#endif // !GA_REAL_HPP_
