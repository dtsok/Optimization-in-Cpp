#ifndef NELDER_MEAD_HPP_
#define NELDER_MEAD_HPP_

#include "../common/common_file.hpp"

class NelderMead {
	private:
	size_t maxIterations;
	double acc;
	double real_val;

	double **simplex;
	double (*function)(const size_t, const double *);
	size_t N;

	const double r_inc = -0.5f;
	const double r_exc = 0.5f;
	const double r_ref = 1.0f;
	const double r_exp = 2.0f;

	size_t iterations = 0;
	int inc_counter = 0;
	int exc_counter = 0;
	int ref_counter = 0;
	int exp_counter = 0;
	int shr_counter = 0;
	int reset_counter = 0;

	double *values;
	double *cm;
	double *x_inc;
	double *x_ref;
	double *x_exc;
	double *x_exp;

	double *best_point;
	double g_minimum;

	void centerMass();
	void generatePoint(double *x_op, const double r_op);
	void shrinkSimplex(double beta);

	int partition(int l, int h);
	void quicksort(int l, int h);

	void checkBeforeFree(double *other);
	void resetSimplex();

	public:
	NelderMead(size_t dimensions, double **points, double (*func)(const size_t, const double *));
	void setParameters(size_t iter, double acc, double real);
	~NelderMead();

	void minimize();
};

#endif // NELDER_MEAD_HPP_