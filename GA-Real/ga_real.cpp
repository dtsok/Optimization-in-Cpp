#include "ga_real.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <limits>


GA_Real::GA_Real(size_t dimensions, double (*func)(const size_t, const double *), double l, double h){
	dim = dimensions;
	function = func;
	this->l=l;
	this->h=h;
}

GA_Real::~GA_Real()
{
	if (P != nullptr) {
		for (size_t i = 0; i < populationSize; i++) {
			delete[] P[i];
		}
		delete[] P;
	}

	if (S != nullptr) {
		for (size_t i = 0; i < populationSize; i++) {
			delete[] S[i];
		}
		delete[] S;
	}
	if (p_values != nullptr) {
		delete[] p_values;
	}

	if (best != nullptr) {
		delete[] best;
	}

	if (worst != nullptr) {
		delete[] worst;
	}

	if (s_values != nullptr) {
		delete[] s_values;
	}
}


void GA_Real::setParameters(size_t iterations, double acc, double real_val, size_t popSize)
{
	maxIterations = iterations;
	this->acc = acc;
	this->real_val = real_val;
	populationSize = popSize;

}

void GA_Real::initialize()
{
	P = new double*[populationSize];
	p_values = new double[populationSize]();

	S = new double*[populationSize];
	s_values = new double[populationSize]();

	best = new double[dim]();
	worst = new double[dim]();

	for (size_t i = 0; i < populationSize; i++) {
		P[i] = new double[dim];
		S[i] = new double[dim];
		for (size_t j = 0; j < dim; j++) {
			P[i][j] = RandomGenerator::generateDouble(l, h);
			S[i][j] = 0;
		}
	}
}

void GA_Real::evaluate(double **pop, double *val, size_t N)
{
	worst_value = std::numeric_limits<double>::min();
	for (size_t i = 0; i < N; i++) {
		val[i] = function(dim, pop[i]);
		if (val[i] < best_value) {
			best_value = val[i];
			// std::copy(pop[i], pop[i] + total_bits, best);
		}
		if (val[i] > worst_value) {
			worst_value = val[i];
		
		}
	}
	// std::cout<<"best "<<best_value<<" worst "<<worst_value<<"\n";
}

void GA_Real::tournament_selection(size_t tourSize)
{
	std::set<size_t> members;
	for (size_t i = 0; i < populationSize; i++) {
		members.clear();
		while (members.size() < tourSize) {
			members.insert(RandomGenerator::generateInt(0, populationSize - 1));
		}
		size_t best_index = 0;
		double best_v = std::numeric_limits<double>::max();
		for (const auto &mem : members) {
			if (p_values[best_index] < best_v) {
				best_v = p_values[best_index];
				best_index = mem;
			}
		}

		std::copy(P[best_index], P[best_index] + dim, S[i]);
	}
}

void GA_Real::nonlinear_ranking(double *l_bounds)
{
	double *fitness = new double[populationSize]();
	double total_fitness = 0;
	for (size_t i = 0; i < populationSize; i++) {
		if (p_values[i]>worst_value) {
			worst_value = p_values[i];
		}
	}
	for (size_t i = 0; i < populationSize; i++) {
		fitness[i] = worst_value - p_values[i];
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

void GA_Real::linear_ranking(double *l_bounds)
{
	double sp = 1.5; // selective pressure
	double *fitness = new double[populationSize]();
	double total_fitness = 0;
	int *sortedInd = new int[populationSize];
	for (size_t i = 0; i < populationSize; i++) {
		sortedInd[i] = i;
	}
	std::sort(sortedInd, sortedInd + populationSize, [&](const int &a, const int &b) { return (p_values[a] > p_values[b]); });

	for (size_t i = 0; i < populationSize; i++) {
		fitness[sortedInd[i]] = 2 - sp + 2 * (sp - 1) * (i) / (populationSize - 1);
		total_fitness += fitness[sortedInd[i]];
	}

	fitness[0] = fitness[0] / total_fitness;
	for (size_t i = 1; i < populationSize; i++) {
		fitness[i] = fitness[i] / total_fitness;
		l_bounds[i] = l_bounds[i - 1] + fitness[i - 1];
	}

	delete[] sortedInd;
	delete[] fitness;
}

void GA_Real::roulette_wheel_selection(bool linear)
{
	double *l_bounds = new double[populationSize]();

	if (linear) {
		linear_ranking(l_bounds);
	}
	else {
		nonlinear_ranking(l_bounds);
	}

	size_t index = 0;
	while (index < populationSize) {
		double r = RandomGenerator::generateDouble(0, 1);
		for (size_t j = 1; j < populationSize; j++) {
			if (r <= l_bounds[j]) {
				std::copy(P[j - 1], P[j - 1] + dim, S[index]);
				index++;
				break;
			}
		}
	}

	// evaluate(selected, populationSize);
	delete[] l_bounds;
}

void GA_Real::crossoverAndUpdate(const size_t p1, const size_t p2, const double delta){
	double *descendant = new double[dim]();
	for (size_t i = 0; i < dim; i++) {
		double r = RandomGenerator::generateDouble(-delta, 1+delta);
		descendant[i] = r*S[p1][i] + (1-r)*S[p2][i];
		if (descendant[i] < l) {
			descendant[i] = l;
		}
		else if (descendant[i] > h) {
			descendant[i] = h;
		}
	}
	if (p_values[p1] < p_values[p2]) {
		std::copy(descendant, descendant+dim, S[p2]);
	}
	else {
		std::copy(descendant, descendant+dim, S[p1]);
	}

	delete [] descendant;
}

void GA_Real::crossover(double crop){
	double delta = 0.15; // used in order to avoid gradual shrinkage
	
	std::vector<size_t> parents;
	for (size_t i = 0; i < populationSize; i++) {
		double r = RandomGenerator::generateDouble(0, 1);
		if (r <= crop) {
			parents.push_back(i);
		}
	}

	if (parents.size() % 2 != 0 && parents.size() == populationSize) {
		parents.pop_back();
	}

	while (parents.size() % 2 != 0 && parents.size() < populationSize) {
		parents.push_back(RandomGenerator::generateInt(0, populationSize - 1));
	}

	for (size_t i = 0; i < parents.size(); i += 2) {
		crossoverAndUpdate(parents[i], parents[i + 1], delta);
	}

}

void GA_Real::mutation(double mutop){
	for (size_t i = 0; i < populationSize; i++) {
		for (size_t j = 0; j < dim; j++) {
			if (RandomGenerator::generateDouble(0, 1)<mutop) {
				double alpha = std::min(l+std::abs(S[i][j]), h - std::abs(S[i][j]));
				double z = RandomGenerator::generateDouble(-alpha, alpha);
				S[i][j] += z;
				if (S[i][j] < l) {
					S[i][j] = l;
				}
				else if (S[i][j] > h) {
					S[i][j] = h;
				}
			}
		}
	}
}

void GA_Real::newPopulation()
{
	// double *selected_values = new double[populationSize]();
	// evaluate(S, s_values, populationSize);

	int *indicies_pop = new int[populationSize]();
	int *indicies_selected = new int[populationSize]();
	for (size_t i = 0; i < populationSize; i++) {
		indicies_pop[i] = indicies_selected[i] = i;
	}

	std::sort(indicies_pop, indicies_pop + populationSize, [&](const int &a, const int &b) { return (p_values[a] < p_values[b]); });

	std::sort(indicies_selected, indicies_selected + populationSize, [&](const int &a, const int &b) { return (s_values[a] < s_values[b]); });

	double **temp = new double *[populationSize];
	size_t current = 0;
	size_t ind_1 = 0;
	size_t ind_2 = 0;
	while (current < populationSize) {
		temp[current] = new double[dim]();
		std::copy(P[indicies_pop[ind_1]], P[indicies_pop[ind_1]] + dim, temp[current]);
		ind_1++;
		current++;
		if (current < populationSize) {
			temp[current] = new double[dim]();
			std::copy(S[indicies_selected[ind_2]], S[indicies_selected[ind_2]] + dim, temp[current]);
			ind_2++;
			current++;
		}
	}

	for (size_t i = 0; i < populationSize; i++) {
		std::copy(temp[i], temp[i] + dim, P[i]);
	}

	for (size_t i = 0; i < populationSize; i++) {
		delete[] temp[i];
	}
	delete[] temp;

	delete[] indicies_pop;
	delete[] indicies_selected;
}

void GA_Real::minimize(bool tour)
{
	initialize();
	evaluate(P, p_values, populationSize);
	double crossover_rate = 0.5;
	double mutation_rate = 0.2;
	//	bool flag = true;
	size_t iter = 0;
	while (best_value > real_val + acc && iter < maxIterations) {
		// double prev_ = best_value;
		if (tour) {
			tournament_selection(populationSize / 2);
		}
		else {
			roulette_wheel_selection(true);
		}
			crossover(crossover_rate);
		mutation(mutation_rate);
		evaluate(S, s_values, populationSize);
		newPopulation();
		evaluate(P, p_values, populationSize);
		std::cout << iter << ": " << best_value << "\n";
		iter++;
	}
	std::cout << iter << ": " << best_value << "\n";
}
