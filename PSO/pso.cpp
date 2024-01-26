#include "pso.hpp"

PSO::PSO(size_t dimensions, double (*func)(const size_t, const double *), double l, double h)
{
	dim = dimensions;
	function = func;
	this->l = l;
	this->h = h;
}

PSO::~PSO()
{
	if (swarm != nullptr) {
		for (size_t i = 0; i < swarmSize; i++) {
			delete[] swarm[i];
		}
		delete[] swarm;
	}

	if (neighborhood != nullptr) {
		for (size_t i = 0; i < swarmSize; i++) {
			delete[] neighborhood[i];
		}
		delete[] neighborhood;
	}

	if (velocity != nullptr) {
		for (size_t i = 0; i < swarmSize; i++) {
			delete[] velocity[i];
		}
		delete[] velocity;
	}

	if (bestPositions != nullptr) {
		for (size_t i = 0; i < swarmSize; i++) {
			delete[] bestPositions[i];
		}
		delete[] bestPositions;
	}

	if (s_values != nullptr) {
		delete[] s_values;
	}

	if (pos_values != nullptr) {
		delete[] pos_values;
	}
}

void PSO::setParameters(size_t iterations, double acc, double real_val, size_t swarmSize)
{
	maxIterations = iterations;
	this->acc = acc;
	this->real_val = real_val;
	this->swarmSize = swarmSize;
}

void PSO::initialize()
{
	swarm = new double *[swarmSize];
	s_values = new double[swarmSize]();
	pos_values = new double[swarmSize];

	velocity = new double *[swarmSize];
	bestPositions = new double *[swarmSize];

	for (size_t i = 0; i < swarmSize; i++) {
		swarm[i] = new double[dim];
		velocity[i] = new double[dim];
		bestPositions[i] = new double[dim];
		pos_values[i] = std::numeric_limits<double>::max();
		for (size_t j = 0; j < dim; j++) {
			swarm[i][j] = RandomGenerator::generateDouble(l, h);
			bestPositions[i][j] = swarm[i][j];
			velocity[i][j] = .5;
		}
	}
}

void PSO::checkBounds(double **array)
{
	for (size_t i = 0; i < swarmSize; i++) {
		for (size_t j = 0; j < dim; j++) {
			if (array[i][j] < l) {
				array[i][j] = l;
			}
			else if (array[i][j] > h) {
				array[i][j] = h;
			}
		}
	}
}

void PSO::checkBoundsVelocity(double **array, double alpha)
{
	double vmax = alpha * (h - l);
	for (size_t i = 0; i < swarmSize; i++) {
		for (size_t j = 0; j < dim; j++) {
			if (array[i][j] < -vmax) {
				array[i][j] = -vmax;
			}
			else if (array[i][j] > vmax) {
				array[i][j] = vmax;
			}
		}
	}
}

void PSO::evaluateAndUpdateBest(double **array, double *store, size_t N)
{
	for (size_t i = 0; i < N; i++) {
		double val = function(dim, array[i]);
		if (val < best_value) {
			best_value = val;
		}
		if (val > worst_value) {
			worst_value = val;
		}

		if (val < pos_values[i]) {
			pos_values[i] = val;
			std::copy(array[i], array[i] + dim, bestPositions[i]);
		}
	}
}
void PSO::initializeNeighborhood(bool topology)
{
	neighborhood = new int *[swarmSize];
	radius = topology ? swarmSize : 3;
	for (size_t i = 0; i < swarmSize; i++) {
		neighborhood[i] = new int[radius];
		if (topology) {
			for (size_t j = 0; j < swarmSize; j++) {
				neighborhood[i][j] = j;
			}
		}
		else {
			if (i == 0) {
				neighborhood[i][0] = swarmSize - 1;
			}
			else {
				neighborhood[i][0] = i - 1;
			}
			neighborhood[i][1] = i;
			neighborhood[i][2] = (i + 1) % swarmSize;
		}
	}
}

void PSO::updateSwarm()
{
	for (size_t i = 0; i < swarmSize; i++) {
		for (size_t j = 0; j < dim; j++) {
			swarm[i][j] += velocity[i][j];
		}
	}
}

size_t PSO::findBest(size_t i)
{
	size_t indx = 0;
	double val = std::numeric_limits<double>::max();
	for (size_t j = 0; j < radius; j++) {
		if (pos_values[neighborhood[i][j]] < val) {
			val = pos_values[neighborhood[i][j]];
			indx = neighborhood[i][j];
		}
	}
	return indx;
}

void PSO::updateVelocity()
{
	for (size_t i = 0; i < swarmSize; i++) {
		for (size_t j = 0; j < dim; j++) {
			double r1 = RandomGenerator::generateDouble(0, 1);
			double r2 = RandomGenerator::generateDouble(0, 1);

			size_t best_from_ng = findBest(i);
			velocity[i][j] =
				xi * (velocity[i][j] + r1 * c1 * (bestPositions[i][j] - swarm[i][j]) + r2 * c2 * (bestPositions[best_from_ng][j] - swarm[i][j]));
		}
	}
}
void PSO::resetSwarm()
{
	for (size_t i = 0; i < swarmSize; i++) {
		for (size_t j = 0; j < dim; j++) {
			swarm[i][j] = RandomGenerator::generateDouble(l, h);
			velocity[i][j] = .8;
		}
	}
	evaluateAndUpdateBest(swarm, s_values, swarmSize);
}

void PSO::minimize(bool topology)
{
	initialize();
	initializeNeighborhood(topology);
	evaluateAndUpdateBest(swarm, s_values, swarmSize);
	size_t iter = 0;
	while (best_value > real_val + acc && iter < maxIterations) {
		updateVelocity();
		checkBoundsVelocity(velocity);
		updateSwarm();
		checkBounds(swarm);
		evaluateAndUpdateBest(swarm, s_values, swarmSize);
		std::cout << iter << ": " << best_value << "\n";
		iter++;
	}
	std::cout << iter << ": " << best_value << "\n";
}
