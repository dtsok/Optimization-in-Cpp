#include <iostream>
#include <random>

std::default_random_engine generator;
std::uniform_real_distribution<double> dist(-2.5, 2.5);

/* With arrays */

double l2_norm(const double *X, const int i, const int j)
{
	double val = 0;
	for (int z = 0; z < 3; z++) {
		val += (X[i + z] - X[j + z]) * (X[i + z] - X[j + z]);
	}

	return std::sqrt(val);
}

double ELJ(const int N, const double *X, const double epsilon = 1, const double sigma = 1)
{
	double sum = 0;
	for (size_t i = 0; i < N - 3; i += 3) {
		double temp = 0;
		for (size_t j = i + 3; j < N; j += 3) {
			double rij = l2_norm(X, i, j);
			double factor = std::pow((sigma / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * epsilon * sum;
}

void centerMass(const int N, double **points, double *cm)
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

int partition(double **points, double *values, int l, int h)
{
	// int R = rand() % (h - l + 1) + l;
	int R = h;
	double pivot = values[R];
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

void quicksort(double **points, double *values, int l, int h)
{
	if (l < h) {

		int pivot_part = partition(points, values, l, h);

		quicksort(points, values, l, pivot_part - 1);
		quicksort(points, values, pivot_part + 1, h);
	}
}

void generate_point(double *x_op, const double *worst, double *cm, const int N, const double r_op)
{
	for (size_t i = 0; i < N; i++) {
		x_op[i] = (1 + r_op) * cm[i] - r_op * worst[i];
	}
}

void shrinkSimplex(const int N, double **points, double *values, double beta = .5)
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

bool checkBeforeFree(double *point2delete, double *x1, double *x2, double *x3, double *x4)
{
	if (point2delete == x1 || point2delete == x2 || point2delete == x3 || point2delete == x4) {
		return false;
	}
	return true;
}

void freePoints(double **points, double *other, const int N)
{
	for (size_t i = 0; i < N + 1; i++) {
		if (points[i] == other) {
			points[i] = nullptr;
			// break;
		}
	}
}

void NelderMead(const int N, double **points)
{
	double *values = new double[N + 1]();
	double *cm = new double[N]();
	double *x_inc = new double[N]();
	double *x_ref = new double[N]();
	double *x_exc = new double[N]();
	double *x_exp = new double[N]();

	for (size_t i = 0; i < N + 1; i++) {
		values[i] = ELJ(N, points[i]);
	}
	quicksort(points, values, 0, N);

	const double r_inc = -0.5;
	const double r_exc = 0.5;
	const double r_ref = 1;
	const double r_exp = 2;

	int inc_counter = 0;
	int exc_counter = 0;
	int ref_counter = 0;
	int exp_counter = 0;

	int iterations = 0;
	const int maxIterations = 500000;
	const double real_val = -9.10385;
	const double acc = 0.001;
	double of_value = values[0];
	double g_minimum = values[0];
	while (of_value > acc + real_val && iterations < maxIterations) {
		centerMass(N, points, cm);
		generate_point(x_ref, points[N], cm, N, r_ref);
		double x_ref_val = ELJ(N, x_ref);
		if (x_ref_val >= values[0] && x_ref_val < values[N - 1]) {
			if (checkBeforeFree(points[N], x_ref, x_exp, x_exc, x_inc)) {
				delete[] points[N];
			}
			points[N] = x_ref;
			values[N] = x_ref_val;
		}
		else if (x_ref_val < values[0]) {
			generate_point(x_exp, points[N], cm, N, r_exp);
			double x_exp_val = ELJ(N, x_exp);
			if (x_exp_val < x_ref_val) {
				if (checkBeforeFree(points[N], x_ref, x_exp, x_exc, x_inc)) {
					delete[] points[N];
				}
				points[N] = x_exp;
				values[N] = x_exp_val;
			}
			else {
				if (checkBeforeFree(points[N], x_ref, x_exp, x_exc, x_inc)) {
					delete[] points[N];
				}
				points[N] = x_ref;
				values[N] = x_ref_val;
			}
		}
		else if (x_ref_val >= values[N - 1] && x_ref_val < values[N]) {
			generate_point(x_exc, points[N], cm, N, r_exc);
			double x_exc_val = ELJ(N, x_exc);
			if (x_exc_val <= x_ref_val) {
				if (checkBeforeFree(points[N], x_ref, x_exp, x_exc, x_inc)) {
					delete[] points[N];
				}
				points[N] = x_exc;
				values[N] = x_exc_val;
			}
			else {
				if (checkBeforeFree(points[N], x_ref, x_exp, x_exc, x_inc)) {
					delete[] points[N];
				}
				points[N] = x_ref;
				values[N] = x_ref_val;
			}
		}
		else if (values[N] <= x_ref_val) {
			generate_point(x_inc, points[N], cm, N, r_inc);
			double x_inc_val = ELJ(N, x_inc);
			if (x_inc_val < values[N]) {
				if (checkBeforeFree(points[N], x_ref, x_exp, x_exc, x_inc)) {
					delete[] points[N];
				}
				points[N] = x_inc;
				values[N] = x_inc_val;
			}
			else {
				shrinkSimplex(N, points, values);
			}
		}

		quicksort(points, values, 0, N);
		of_value = values[0];
		if (of_value < g_minimum) {
			g_minimum = of_value;
		}

		std::cout << iterations << ": " << of_value << " global: " << g_minimum << "\n";
		iterations++;
	}

	delete[] values;
	delete[] cm;

	freePoints(points, x_inc, N);
	freePoints(points, x_ref, N);
	freePoints(points, x_exc, N);
	freePoints(points, x_exp, N);
	delete[] x_inc;
	// x_inc = nullptr;
	delete[] x_ref;
	// x_ref = nullptr;
	delete[] x_exc;
	// x_exc = nullptr;
	delete[] x_exp;
	// x_exp = nullptr;
	// return ;
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
	double **p = new double *[dim + 1]; //{0.7430002202, 0.2647603899, -0.0468575389, -0.7430002647, -0.2647604843, 0.0468569750, 0.1977276118,
										//-0.4447220146, 0.6224700350,
										// -0.1977281310, 0.4447221826, -0.6224697723, -0.1822009635, 0.5970484122, 0.4844363476, 0.1822015272,
										// -0.5970484858, -0.4844360463};

	for (size_t t = 0; t < dim + 1; t++) {
		p[t] = new double[dim];
		for (size_t z = 0; z < dim; z++) {
			p[t][z] = dist(generator);
			// std::cout<<p[t][z]<<" ";
		}
		// std::cout<<"\n";
	}

	NelderMead(dim, p);

	// delete[] points;
	for (size_t t = 0; t < dim + 1; t++) {
		if (p[t] != nullptr) {
			delete[] p[t];
			p[t] = nullptr;
		}
	}
	delete[] p;
	return 0;
}

// for x = 5
/*double x[] = {-0.2604720088	,	0.7363147287	,	0.4727061929	,
0.260471655	,	-0.7363150782	,	-0.4727063011	,
-0.4144908003	,	-0.3652598516	,	0.340555962	,
-0.1944131041	,	0.2843471802	,	-0.5500413671	,
0.6089042582	,	0.0809130209	,	0.2094855133
};
std::cout<<"ANS "<<ELJ(5, x)<<"\n";*/