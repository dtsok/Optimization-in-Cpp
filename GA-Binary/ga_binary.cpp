#include "ga_binary.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <set>
#include <vector>

GA_Binary::~GA_Binary()
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

GA_Binary::GA_Binary(size_t dimensions, double (*func)(const size_t, const double *), double l, double h)
{
	dim = dimensions;
	function = func;
	this->l = l;
	this->h = h;
}

void GA_Binary::setParameters(size_t iterations, double acc, double real_val, size_t popSize)
{
	maxIterations = iterations;
	this->acc = acc;
	this->real_val = real_val;
	populationSize = popSize;

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

void GA_Binary::initialize()
{
	P = new int *[populationSize];
	p_values = new double[populationSize]();

	S = new int *[populationSize];
	s_values = new double[populationSize]();

	best = new int[total_bits]();
	worst = new int[total_bits]();

	for (size_t i = 0; i < populationSize; i++) {
		P[i] = new int[total_bits];
		S[i] = new int[total_bits];
		for (size_t j = 0; j < total_bits; j++) {
			P[i][j] = RandomGenerator::generateInt(0, 1);
			S[i][j] = 0;
		}
	}
}

void GA_Binary::evaluate(int **pop, double *val, size_t N)
{
	double *storeRealVector = new double[dim]();
	for (size_t i = 0; i < N; i++) {
		size_t index = 0;
		for (size_t j = 0; j < total_bits; j += bits_per_num) {
			storeRealVector[index++] = decode(pop[i], j);
		}
		val[i] = function(dim, storeRealVector);
		if (val[i] < best_value) {
			best_value = val[i];
			// std::copy(pop[i], pop[i] + total_bits, best);
		}
		else if (val[i] > worst_value) {
			worst_value = val[i];
			// std::copy(pop[i], pop[i] + total_bits, worst);
		}
		// std::cout << "values = " << val[i] << "\n";
	}
	// std::cout<<"best "<<best_value<<" worst "<<worst_value<<"\n";
	delete[] storeRealVector;
}

void GA_Binary::tournament_selection(size_t tourSize)
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

		std::copy(P[best_index], P[best_index] + total_bits, S[i]);
	}
}

void GA_Binary::nonlinear_ranking(double *l_bounds)
{
	double *fitness = new double[populationSize]();
	double total_fitness = 0;
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

void GA_Binary::linear_ranking(double *l_bounds)
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

void GA_Binary::roulette_wheel_selection(bool linear)
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
				std::copy(P[j - 1], P[j - 1] + total_bits, S[index]);
				index++;
				break;
			}
		}
	}

	// evaluate(selected, populationSize);
	delete[] l_bounds;
}

void GA_Binary::crossoverAndUpdate(const size_t p1, const size_t p2, const std::set<size_t> &ind)
{
	bool flag = true;
	size_t current = 0; // starting index
	for (const size_t &index : ind) {
		while (current < index) {
			if (flag) {
				std::swap(S[p1][current], S[p2][current]);
			}
			current++;
		}
		flag = flag ? false : true;
	}

	// swap rest elements if needed
	if ((current < total_bits - 1) && flag) {
		for (size_t i = current; i < total_bits; i++) {
			std::swap(S[p1][i], S[p2][i]);
		}
	}
}

void GA_Binary::crossover(double crop)
{
	// number of crossover points and corresponding places
	size_t k = dim;
	// size_t k = 10;
	std::set<size_t> indicies;
	for (size_t i = 0; i < k; i++) {
		indicies.insert(RandomGenerator::generateInt(1, total_bits - 2));
	}

	std::vector<size_t> parents;
	for (size_t i = 0; i < populationSize; i++) {
		double r = RandomGenerator::generateDouble(0, 1);
		if (r <= crop) {
			parents.push_back(i);
		}
	}

	if (parents.size() % 2 != 0) {
		parents.push_back(RandomGenerator::generateInt(0, populationSize - 1));
	}

	for (size_t i = 0; i < parents.size(); i += 2) {
		crossoverAndUpdate(parents[i], parents[i + 1], indicies);
	}
}

void GA_Binary::mutation(double mutop)
{
	for (size_t i = 0; i < populationSize; i++) {
		for (size_t j = 0; j < total_bits; j++) {
			if (RandomGenerator::generateDouble(0, 1) < mutop) {
				if (S[i][j]) {
					S[i][j] = 0;
				}
				else {
					S[i][j] = 1;
				}
			}
		}
	}
}

void GA_Binary::newPopulation()
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

	int **temp = new int *[populationSize];
	size_t current = 0;
	size_t ind_1 = 0;
	size_t ind_2 = 0;
	while (current < populationSize) {
		temp[current] = new int[total_bits]();
		if (p_values[indicies_pop[ind_1]] < s_values[indicies_selected[ind_2]] && ind_1 < total_bits) {
			std::copy(P[indicies_pop[ind_1]], P[indicies_pop[ind_1]] + total_bits, temp[current]);
			ind_1++;
		}
		else if (p_values[indicies_pop[ind_1]] > s_values[indicies_selected[ind_2]] && ind_2 < total_bits) {
			std::copy(S[indicies_selected[ind_2]], S[indicies_selected[ind_2]] + total_bits, temp[current]);
			ind_2++;
		}
		current++;
	}

	for (size_t i = 0; i < populationSize; i++) {
		std::copy(temp[i], temp[i] + total_bits, P[i]);
	}

	for (size_t i = 0; i < populationSize; i++) {
		delete[] temp[i];
	}
	delete[] temp;

	delete[] indicies_pop;
	delete[] indicies_selected;
}

void GA_Binary::minimize(bool tour)
{
	initialize();
	evaluate(P, p_values, populationSize);
	double crossover_rate = 0.2;
	double mutation_rate = 0.1;
	//	bool flag = true;
	size_t iter = 0;
	while (best_value > real_val + acc && iter < maxIterations) {
		// double prev_ = best_value;
		if (tour) {
			tournament_selection(populationSize / 2);
		}
		else {
			roulette_wheel_selection(false);
		}
		/*if (best_value + std::abs(real_val) < 0.05 && flag) {
			crossover_rate = 0.25;
			mutation_rate = 0.35;
			flag = false;
		}*/
		crossover(crossover_rate);
		mutation(mutation_rate);
		evaluate(S, s_values, populationSize);
		newPopulation();
		evaluate(P, p_values, populationSize);
		// if (best_value < prev_) {
		std::cout << iter << ": " << best_value << "\n";
		//	prev_ = best_value;
		//}
		iter++;
	}
	std::cout << iter << ": " << best_value << "\n";
}
