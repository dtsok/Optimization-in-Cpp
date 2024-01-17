#include <algorithm>
#include <iostream>
#include <random>

std::default_random_engine generator(time(0));
std::uniform_real_distribution<float> dist(-2.5, 2.5);

const int maxIterations = 10000;
const float real_val = -12.712f;
const float acc = 0.001f;
int reset_counter = 0;

float l2_norm(const float *X, const int i, const int j)
{
	float val = 0;
	for (int z = 0; z < 3; z++) {
		val += (X[i + z] - X[j + z]) * (X[i + z] - X[j + z]);
	}

	return std::sqrt(val);
}

float l2_norm_gen(float *x, float *y, const int N)
{
	float val = 0;
	for (size_t i = 0; i < N; i++) {
		val += (x[i] - y[i]) * (x[i] - y[i]);
	}
	return std::sqrt(val);
}

float ELJ(const int N, const float *X, const float epsilon = 1, const float sigma = 1)
{
	float sum = 0;
	for (size_t i = 0; i < N - 3; i += 3) {
		float temp = 0;
		for (size_t j = i + 3; j < N; j += 3) {
			float rij = l2_norm(X, i, j);
			float factor = std::pow((sigma / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * epsilon * sum;
}

void centerMass(const int N, float **points, float *cm)
{
	for (size_t i = 0; i < N; i++) {
		cm[i] = 0;
	}
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			cm[j] += points[i][j] / N;
		}
	}
}

int partition(float **points, float *values, int l, int h)
{
	int R = h;
	float pivot = values[R];
	int temp_right = l - 1;
	for (int i = l; i <= h; i++) {
		if (values[i] < pivot) {
			temp_right++;
			std::swap(values[i], values[temp_right]);
			std::swap(points[temp_right], points[i]);
		}
	}

	std::swap(values[temp_right + 1], values[h]);
	std::swap(points[temp_right + 1], points[h]);

	return temp_right + 1;
}

void quicksort(float **points, float *values, int l, int h)
{
	if (l < h) {

		int pivot_part = partition(points, values, l, h);

		quicksort(points, values, l, pivot_part - 1);
		quicksort(points, values, pivot_part + 1, h);
	}
}

void generate_point(float *x_op, const float *worst, float *cm, const int N, const float r_op)
{

	for (size_t i = 0; i < N; i++) {
		x_op[i] = (1 + r_op) * cm[i] - r_op * worst[i];
	}
}

void shrinkSimplex(const int N, float **points, float *values, float beta = .5)
{
	for (size_t i = 1; i < N + 1; i++) {
		for (size_t j = 0; j < N; j++) {
			points[i][j] = points[0][j] - (points[i][j] - points[0][j]) * beta;
		}
	}

	for (size_t i = 1; i < N + 1; i++) {
		values[i] = ELJ(N, points[i]);
	}
}

void checkBeforeFree(float **points, float *other, const int N)
{
	for (size_t i = 0; i < N + 1; i++) {
		if (points[i] == other) {
			points[i] = nullptr;
			break;
		}
	}
}

void resetSimplex(float **points, const int N)
{
	if (l2_norm_gen(points[0], points[N], N) < 0.01) {
		reset_counter++;
		for (size_t t = 0; t < N + 1; t++) {
			for (size_t z = 0; z < N; z++) {
				points[t][z] = dist(generator);
			}
		}
	}
}

void NelderMead(const int N, float **points)
{
	float *values = new float[N + 1]();
	float *cm = new float[N]();
	float *x_inc = new float[N]();
	float *x_ref = new float[N]();
	float *x_exc = new float[N]();
	float *x_exp = new float[N]();
	float *best_point = new float[N]();

	for (size_t i = 0; i < N + 1; i++) {
		values[i] = ELJ(N, points[i]);
	}
	quicksort(points, values, 0, N);

	const float r_inc = -0.5f;
	const float r_exc = 0.5f;
	const float r_ref = 1.0f;
	const float r_exp = 2.0f;

	int inc_counter = 0;
	int exc_counter = 0;
	int ref_counter = 0;
	int exp_counter = 0;
	int shr_counter = 0;

	int iterations = 0;

	float of_value = values[0];
	float g_minimum = values[0];
	while (of_value > acc + real_val && iterations < maxIterations) {
		resetSimplex(points, N);
		centerMass(N, points, cm);
		generate_point(x_ref, points[N], cm, N, r_ref);
		float x_ref_val = ELJ(N, x_ref);
		if (x_ref_val >= values[0] && x_ref_val < values[N - 1]) {
			std::copy(x_ref, x_ref + N, points[N]);
			values[N] = x_ref_val;
			ref_counter++;
		}
		else if (x_ref_val < values[0]) {
			generate_point(x_exp, points[N], cm, N, r_exp);
			float x_exp_val = ELJ(N, x_exp);
			if (x_exp_val < x_ref_val) {
				std::copy(x_exp, x_exp + N, points[N]);
				values[N] = x_exp_val;
				exp_counter++;
			}
			else {
				std::copy(x_ref, x_ref + N, points[N]);
				values[N] = x_ref_val;
				ref_counter++;
			}
		}
		else if (x_ref_val >= values[N - 1] && x_ref_val < values[N]) {
			generate_point(x_exc, points[N], cm, N, r_exc);
			float x_exc_val = ELJ(N, x_exc);
			if (x_exc_val <= x_ref_val) {
				std::copy(x_exc, x_exc + N, points[N]);
				values[N] = x_exc_val;
				exc_counter++;
			}
			else {
				std::copy(x_ref, x_ref + N, points[N]);
				values[N] = x_ref_val;
				ref_counter++;
			}
		}
		else if (values[N] <= x_ref_val) {
			generate_point(x_inc, points[N], cm, N, r_inc);
			float x_inc_val = ELJ(N, x_inc);
			if (x_inc_val < values[N]) {
				std::copy(x_inc, x_inc + N, points[N]);
				values[N] = x_inc_val;
				inc_counter++;
			}
			else {
				shrinkSimplex(N, points, values);
				shr_counter++;
			}
		}

		quicksort(points, values, 0, N);
		of_value = values[0];
		if (of_value < g_minimum) {
			g_minimum = of_value;
			std::copy(points[0], points[0] + N, best_point);
		}

		std::cout << iterations << ": " << of_value << " global: " << g_minimum << "\n";
		iterations++;
	}
	std::cout << "counters:\nreflections " << ref_counter << "\nexpansions " << exp_counter << "\ncontractions " << exc_counter
			  << "\ninner contractions " << inc_counter << "\nresets " << reset_counter << "\n";
	// for (size_t i = 0; i < N; i++) {
	// 	std::cout << best_point[i] << ", ";
	// }
	// std::cout << "\n";
	// std::cout << ELJ(N, best_point) << "\n";

	delete[] values;
	delete[] cm;

	checkBeforeFree(points, x_inc, N);
	checkBeforeFree(points, x_ref, N);
	checkBeforeFree(points, x_exc, N);
	checkBeforeFree(points, x_exp, N);

	delete[] x_inc;
	delete[] x_ref;
	delete[] x_exc;
	delete[] x_exp;
	delete[] best_point;
}

int main(int argc, char const *argv[])
{
	int N = 1;
	if (argc > 1) {
		N = std::stoi(argv[1]);
	}
	else {
		return 0;
	}

	int dim = 3 * N;
	float **p = new float *[dim + 1];
	for (size_t t = 0; t < dim + 1; t++) {
		p[t] = new float[dim];
		for (size_t z = 0; z < dim; z++) {
			p[t][z] = dist(generator);
		}
	}

	NelderMead(dim, p);

	for (size_t t = 0; t < dim + 1; t++) {
		if (p[t] != nullptr) {
			delete[] p[t];
		}
	}
	delete[] p;
	return 0;
}
