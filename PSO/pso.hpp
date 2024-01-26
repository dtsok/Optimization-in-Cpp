#ifndef PSO_HPP_
#define PSO_HPP_

#include "../common/common_file.hpp"

class PSO {
	public:
	PSO(size_t dimensions, double (*func)(const size_t, const double *), double l = -2.5, double h = 2.5);
	PSO();
	~PSO();
	void setParameters(size_t iterations, double acc, double real_val, size_t popSize);
	void minimize(bool topology);

	public:
	size_t dim = 0;									  // dimension
	double (*function)(const size_t, const double *); // objective function
	double l = 0;									  // low for search-space [l, h]
	double h = 0;									  // high for search-space [l, h]

	size_t maxIterations = 0;
	double acc = 0;
	double real_val = 0; // (global) known minimum

	/* ~~~ default values for next 3 parameters ~~~ */
	const double c1 = 2.05;	 // cognitive constant
	const double c2 = 2.05;	 // social constant
	const double xi = 0.729; // constriction coefficient
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

	size_t swarmSize = 0;
	double **swarm = nullptr;
	double *s_values = nullptr;

	int **neighborhood = nullptr;

	double **velocity = nullptr;
	double **bestPositions = nullptr;
	double *pos_values = nullptr;

	double best_value = std::numeric_limits<double>::max();
	double worst_value = std::numeric_limits<double>::min();

	void initialize();
	size_t radius = 0;
	void initializeNeighborhood(bool topology);

	void checkBounds(double **array);
	void checkBoundsVelocity(double **array, double alpha = 0.1);

	void evaluateAndUpdateBest(double **pop, double *val, size_t N);

	size_t findBest(size_t i);
	void updateVelocity();
	void updateSwarm();
	void updateBestPositions();
	void resetSwarm();
};

#endif // PSO_HPP_
