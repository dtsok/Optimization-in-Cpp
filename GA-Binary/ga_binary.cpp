#include "ga_binary.hpp"
#include <algorithm>

GA_Binary::~GA_Binary()
{
	if (P != nullptr) {
		for (size_t i = 0; i < populationSize; i++) {
			delete[] P[i];
		}
		delete[] P;
		P = nullptr;
	}
	if (best != nullptr) {
		delete[] best;
	}
	if (worst != nullptr) {
		delete[] worst;
	}
	if (values != nullptr) {
		delete[] values;
	}
}

GA_Binary::GA_Binary(size_t dimensions, double (*func)(const size_t, const double *), double l, double h)
{
	dim = dimensions;
	function = func;
	this->l = l;
	this->h = h;
}

void GA_Binary::setParameters(size_t iterations, double acc, double real_val)
{
	maxIterations = iterations;
	this->acc = acc;
	this->real_val = real_val;

	// total numbers that can be represented in the range (l,h) with step=acc
	size_t numbers_in_lh = 1 + (h - l) / acc;
	// number of bits needed for the representation of "number_in_lh" numbers
	size_t high = 2;
	while (numbers_in_lh > high) {
		total_bits++;
		high *= 2;
	}
	total_bits++;
	bits_per_num = total_bits;
	// total_bits = dimensions*bits for 1 number
	total_bits *= dim;
}

double GA_Binary::decode(int *point, size_t offset) // convert binary to real
{
	double z = 0;
	for (size_t i = 0; i < bits_per_num; i++) {
		z += point[i + offset] * std::pow(2, i);
	}
	return l + z * (h - l) / (std::pow(2, bits_per_num) - 1);
}

void GA_Binary::initialize(size_t populationSize)
{
	this->populationSize = populationSize;
	P = new int *[populationSize];
	values = new double[populationSize]();
	best = new int[total_bits]();
	worst = new int[total_bits]();
	for (size_t i = 0; i < populationSize; i++) {
		P[i] = new int[total_bits];
		for (size_t j = 0; j < total_bits; j++) {
			P[i][j] = RandomGenerator::generateInt(0, 1);
		}
	}
}

void GA_Binary::evaluate(int **pop, size_t N)
{
	double *storeRealVector = new double[dim]();
	for (size_t i = 0; i < N; i++) {
		int index = 0;
		for (size_t j = 0; j < total_bits; j += bits_per_num) {
			storeRealVector[index++] = decode(pop[i], j);
		}
		values[i] = function(dim, storeRealVector);
		if (values[i] < best_value) {
			best_index = i;
			best_value = values[i];
			std::copy(pop[i], pop[i] + total_bits, best);
		}
		else if (values[i] > worst_value) {
			worst_index = i;
			worst_value = values[i];
			std::copy(pop[i], pop[i] + total_bits, worst);
		}
		std::cout << "values = " << values[i] << "\n";
	}
	// std::cout<<"best "<<best_value<<" worst "<<worst_value<<"\n";
	delete[] storeRealVector;
}

void GA_Binary::nonlinear_ranking(double *l_bounds)
{
	double *fitness = new double[populationSize]();
	double total_fitness = 0;
	for (size_t i = 0; i < populationSize; i++) {
		fitness[i] = worst_value - values[i];
		total_fitness += fitness[i];
	}

	// convert fitness into probability and store their lower bounds
	fitness[0] = fitness[0] / total_fitness;
	for (size_t i = 1; i < populationSize; i++) {
		fitness[i] = fitness[i] / total_fitness;
		l_bounds[i] = l_bounds[i - 1] + fitness[i - 1];
	}

	delete[] fitness;
}

void GA_Binary::linear_ranking(double *l_bounds)
{
	double *fitness = new double[populationSize]();
	// double total_fitness = 0;
	// TODO: Complete linear_ranking algorithm
	delete[] fitness;
}

void GA_Binary::selection()
{
	int **selected = new int *[populationSize];
	double *l_bounds = new double[populationSize]();

	nonlinear_ranking(l_bounds);

	size_t index = 0;
	while (index < populationSize) {
		double r = RandomGenerator::generateDouble(0, 1);
		for (size_t j = 0; j < populationSize; j++) {
			if (r <= l_bounds[j]) {
				selected[index] = new int[total_bits]();
				std::copy(P[j - 1], P[j - 1] + total_bits, selected[index]);
				index++;
				break;
			}
		}
	}

	std::cout << "New\n";
	evaluate(selected, populationSize);
	for (size_t i = 0; i < populationSize; i++) {
		delete[] selected[i];
	}
	delete[] selected;
	delete[] l_bounds;
}
