#include "nelder_mead.hpp"

NelderMead::NelderMead(const size_t dimensions, double **points, double (*func)(const size_t, const double *))
{
	simplex = points;
	N = dimensions;
	function = func;
}

NelderMead::~NelderMead()
{
	delete[] values;
	delete[] cm;
	delete[] x_inc;
	delete[] x_ref;
	delete[] x_exc;
	delete[] x_exp;
	delete[] best_point;
}

void NelderMead::setParameters(size_t iter, double accuracy, double real)
{
	maxIterations = iter;
	acc = accuracy;
	real_val = real;
}

void NelderMead::centerMass()
{
	for (size_t i = 0; i < N; i++) {
		cm[i] = 0;
	}
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			cm[j] += simplex[i][j] / N;
		}
	}
}

int NelderMead::partition(int l, int h)
{
	int R = h;
	double pivot = values[R];
	int temp_right = l - 1;
	for (int i = l; i <= h; i++) {
		if (values[i] < pivot) {
			temp_right++;
			std::swap(values[i], values[temp_right]);
			std::swap(simplex[temp_right], simplex[i]);
		}
	}

	std::swap(values[temp_right + 1], values[h]);
	std::swap(simplex[temp_right + 1], simplex[h]);

	return temp_right + 1;
}

void NelderMead::quicksort(int l, int h)
{
	if (l < h) {

		int pivot_part = partition(l, h);

		quicksort(l, pivot_part - 1);
		quicksort(pivot_part + 1, h);
	}
}

void NelderMead::generatePoint(double *x_op, const double r_op)
{

	for (size_t i = 0; i < N; i++) {
		x_op[i] = (1 + r_op) * cm[i] - r_op * simplex[N][i];
	}
}

void NelderMead::shrinkSimplex(double beta = .5f)
{
	for (size_t i = 1; i < N + 1; i++) {
		for (size_t j = 0; j < N; j++) {
			simplex[i][j] = simplex[0][j] - (simplex[i][j] - simplex[0][j]) * beta;
		}
	}

	for (size_t i = 1; i < N + 1; i++) {
		// values[i] = ELJ(N, simplex[i]);
		values[i] = function(N, simplex[i]);
	}
}

void NelderMead::checkBeforeFree(double *other)
{
	for (size_t i = 0; i < N + 1; i++) {
		if (simplex[i] == other) {
			simplex[i] = nullptr;
			break;
		}
	}
}

void NelderMead::resetSimplex()
{
	if (l2_norm_gen(simplex[0], simplex[N], N) < 0.01) {
		reset_counter++;
		for (size_t t = 0; t < N + 1; t++) {
			for (size_t z = 0; z < N; z++) {
				simplex[t][z] = RandomGenerator::generateDouble(-2.5, 2.5);
			}
		}
	}
}

void NelderMead::minimize()
{
	values = new double[N + 1]();
	cm = new double[N]();
	x_inc = new double[N]();
	x_ref = new double[N]();
	x_exc = new double[N]();
	x_exp = new double[N]();
	best_point = new double[N]();

	for (size_t i = 0; i < N + 1; i++) {
		// values[i] = ELJ(N, simplex[i]);
		values[i] = function(N, simplex[i]);
	}
	quicksort(0, N);

	double of_value = values[0];
	g_minimum = values[0];
	while (of_value > acc + real_val && iterations < maxIterations) {
		resetSimplex();
		centerMass();
		generatePoint(x_ref, r_ref);
		// double x_ref_val = ELJ(N, x_ref);
		double x_ref_val = function(N, x_ref);
		if (x_ref_val >= values[0] && x_ref_val < values[N - 1]) {
			std::copy(x_ref, x_ref + N, simplex[N]);
			values[N] = x_ref_val;
			ref_counter++;
		}
		else if (x_ref_val < values[0]) {
			generatePoint(x_exp, r_exp);
			// double x_exp_val = ELJ(N, x_exp);
			double x_exp_val = function(N, x_exp);
			if (x_exp_val < x_ref_val) {
				std::copy(x_exp, x_exp + N, simplex[N]);
				values[N] = x_exp_val;
				exp_counter++;
			}
			else {
				std::copy(x_ref, x_ref + N, simplex[N]);
				values[N] = x_ref_val;
				ref_counter++;
			}
		}
		else if (x_ref_val >= values[N - 1] && x_ref_val < values[N]) {
			generatePoint(x_exc, r_exc);
			// double x_exc_val = ELJ(N, x_exc);
			double x_exc_val = function(N, x_exc);
			if (x_exc_val <= x_ref_val) {
				std::copy(x_exc, x_exc + N, simplex[N]);
				values[N] = x_exc_val;
				exc_counter++;
			}
			else {
				std::copy(x_ref, x_ref + N, simplex[N]);
				values[N] = x_ref_val;
				ref_counter++;
			}
		}
		else if (values[N] <= x_ref_val) {
			generatePoint(x_inc, r_inc);
			// double x_inc_val = ELJ(N, x_inc);
			double x_inc_val = function(N, x_inc);
			if (x_inc_val < values[N]) {
				std::copy(x_inc, x_inc + N, simplex[N]);
				values[N] = x_inc_val;
				inc_counter++;
			}
			else {
				shrinkSimplex();
				shr_counter++;
			}
		}

		quicksort(0, N);
		of_value = values[0];
		if (of_value < g_minimum) {
			g_minimum = of_value;
			std::copy(simplex[0], simplex[0] + N, best_point);
		}

		std::cout << iterations << ": " << of_value << " global: " << g_minimum << "\n";
		iterations++;
	}
	std::cout << "~~~~~~~~~~~~~~~~~~~\ncounters:\nreflections " << ref_counter << "\nexpansions " << exp_counter << "\ncontractions "
			  << exc_counter << "\ninner contractions " << inc_counter << "\nresets " << reset_counter << "\n\nminimum found: " << g_minimum
			  << "\nglobal minimun known: " << real_val << "\n";
}