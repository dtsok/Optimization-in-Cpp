#ifndef NELDER_MEAD_HPP_
#define NELDER_MEAD_HPP_

#include "../common_file.hpp"

class NelderMead {
	private:
	float **simplex;
	float (*function)(const int, const float *);
	int N;

	const float r_inc = -0.5f;
	const float r_exc = 0.5f;
	const float r_ref = 1.0f;
	const float r_exp = 2.0f;

	int iterations = 0;
	int inc_counter = 0;
	int exc_counter = 0;
	int ref_counter = 0;
	int exp_counter = 0;
	int shr_counter = 0;
	int reset_counter = 0;

	float *values;
	float *cm;
	float *x_inc;
	float *x_ref;
	float *x_exc;
	float *x_exp;

	float *best_point;
	float g_minimum;

	void centerMass();
	void generatePoint(float *x_op, const float r_op);
	void shrinkSimplex(float beta);

	int partition(int l, int h);
	void quicksort(int l, int h);

	void checkBeforeFree(float *other);
	void resetSimplex();

	public:
	  
    NelderMead( int dimensions, float **points, float (*func)(const int , const float *));
	~NelderMead();

    void minimize();
};

#endif