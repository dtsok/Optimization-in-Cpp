#include <chrono>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <memory>
#include <random>
using namespace std::chrono;

/*double ELJ(int N, const Eigen::MatrixXd &X, double epsilon = 1, double sigma = 1)
{
	double sum = 0;
	for (size_t i = 0; i < N - 1; i++) {
		double temp = 0;
		Eigen::VectorXd vector_1 = X(i, Eigen::all);
		for (size_t j = i + 1; j < N; j++) {
			Eigen::VectorXd vector_2 = X(j, Eigen::all);
			double rij = (vector_1 - vector_2).norm();
			double factor = std::pow((sigma / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * epsilon * sum;
}

inline void swap(Eigen::MatrixXd *&T, const int i, const int j)
{
	Eigen::MatrixXd temp = T[i];
	T[i] = T[j];
	T[j] = temp;
}

std::unique_ptr<Eigen::MatrixXd> centerMass(int N, Eigen::MatrixXd *&points)
{
	std::unique_ptr<Eigen::MatrixXd> cm(new Eigen::MatrixXd(N, 3));
	for (int row = 0; row < N; row++) {
		for (int col = 0; col < 3; col++) {
			(*cm)(row, col) = 0;
		}
	}

	for (int q = 0; q < N; q++) {

		for (int row = 0; row < N; row++) {
			for (int col = 0; col < 3; col++) {
				(*cm)(row, col) += points[q](row, col);
			}
		}
	}
	(*cm) = (*cm) / N;
	return cm;
}

void initialization(int N, Eigen::MatrixXd *&points, double *values)
{
	double min = ELJ(N, points[0]), max_1 = min, max_2 = min;
	int min_index = 0, max_1_index = 0, max_2_index = 0;

	values[0] = min;
	for (size_t i = 1; i < N + 1; i++) {
		values[i] = ELJ(N, points[i]);
	}
	for (int i = 1; i < N + 1; i++) {
		double temp_ = values[i];
		if (temp_ < min) {
			min = temp_;
			min_index = i;
		}
		else if (temp_ >= max_1) {
			max_1 = temp_;
			max_1_index = i;
		}
	}
	swap(points, 0, min_index);
	swap(values, 0, min_index);
	swap(points, N, max_1_index);
	swap(values, N, max_1_index);
	for (int i = 1; i < N; i++) {
		double temp_ = values[i];
		if (temp_ == max_1) {
			max_2 = max_1;
			max_2_index = i;
			break;
		}
		else if (temp_ >= max_2) {
			max_2 = temp_;
			max_2_index = i;
		}
	}

	swap(points, N - 1, max_2_index);
	swap(values, N - 1, max_2_index);
	// delete[] values;
}

int partition(Eigen::MatrixXd *&points, double *a, int l, int h)
{
	// int R = rand() % (h - l + 1) + l;
	int R = h;
	swap(a, h, R);
	swap(points, h, R);
	int pivot = a[h];
	int temp_right = l - 1;
	for (int i = l; i <= h; i++) {
		if (a[i] < pivot) {
			temp_right++;
			swap(a, temp_right, i);
			swap(points, temp_right, i);
		}
	}

	swap(a, temp_right + 1, h);
	swap(points, temp_right + 1, h);

	return temp_right + 1;
}

void quicksort(Eigen::MatrixXd *&points, double *a, int l, int h)
{
	if (l < h) {
		int pivot_part = partition(points, a, l, h);

		// print(a, n);
		// if (pivot_part > 1) {
		quicksort(points, a, l, pivot_part - 1);
		// }
		// if (pivot_part < h - 1) {
		quicksort(points, a, pivot_part + 1, h);
		// }
	}
}

void shrinkSimplex(int N, Eigen::MatrixXd *&points, double *values)
{
	for (size_t i = 1; i < N + 1; i++) {
		points[i] = points[0] - (points[i] - points[0]) * .5;
		// values[i] = ELJ(N, points[i]);
	}
	// initialization(N, points, values);
	quicksort(points, values, 0, N);
}

void NelderMead(int N, Eigen::MatrixXd *&points)
{
	double *values = new double[N + 1];
	for (size_t i = 0; i < N + 1; i++) {
		values[i] = ELJ(N, points[i]);
	}

	// initialization(N, points, values);
	quicksort(points, values, 0, N);

	double r_inc = -.5;
	double r_exc = .5;
	double r_ref = 1;
	double r_exp = 2;

	int iterations = 0;
	int maxIterations = 100000;
	double real_val = -12.712062;
	double acc = 0.001;
	double of_value = ELJ(N, points[0]);
	while (of_value > acc + real_val && iterations < maxIterations) {
		std::unique_ptr<Eigen::MatrixXd> cm = centerMass(N, points);
		Eigen::MatrixXd x_ref = (1 + r_ref) * (*cm) - r_ref * points[N];
		double x_ref_val = ELJ(N, x_ref);
		if (x_ref_val >= values[0] && x_ref_val < values[N - 1]) {
			points[N] = x_ref;
			values[N] = x_ref_val;
		}
		else if (x_ref_val < values[0]) {
			Eigen::MatrixXd x_exp = (1 + r_exp) * (*cm) - r_exp * points[N];
			double x_exp_val = ELJ(N, x_exp);
			if (x_exp_val < x_ref_val) {
				points[N] = x_exp;
				values[N] = x_exp_val;
			}
			else {
				points[N] = x_ref;
				values[N] = x_ref_val;
			}
		}
		else if (x_ref_val >= values[N - 1] && x_ref_val < values[N]) {
			Eigen::MatrixXd x_exc = (1 + r_exc) * (*cm) - r_exc * points[N];
			double x_exc_val = ELJ(N, x_exc);
			if (x_exc_val <= x_ref_val) {
				points[N] = x_exc;
				values[N] = x_exc_val;
			}
		}
		else if (values[N] <= x_ref_val) {
			Eigen::MatrixXd x_inc = (1 + r_inc) * (*cm) - r_inc * points[N];
			double x_inc_val = ELJ(N, x_inc);
			if (x_inc_val < values[N]) {
				points[N] = x_inc;
				values[N] = x_inc_val;
			}
		}
		else {
			shrinkSimplex(N, points, values);
		}
		quicksort(points, values, 0, N);
		of_value = values[0];
		std::cout << iterations << ": " << of_value << "\n";
		iterations++;
	}

	delete[] values;
}
*/

/* With arrays */

double l2_norm(const int N, const double *X, const int i, const int j)
{
	if (std::abs(i - j) < 3) {
		std::cout << "Error @ l2_norm\n";
		exit(-1);
	}
	double val = 0;
	for (int z = 0; z < 3; z++) {
		val += (X[i + z] - X[j + z]) * (X[i + z] - X[j + z]);
	}

	return std::sqrt(val);
}

double ELJ(const int N, const double *X, const double epsilon = 1, const double sigma = 1)
{
	double sum = 0;
	for (size_t i = 0; i < 3 * N - 3; i += 3) {
		double temp = 0;
		for (size_t j = i + 3; j < 3 * N; j += 3) {
			double rij = l2_norm(N, X, i, j);
			double factor = std::pow((sigma / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * epsilon * sum;
}

inline void swap(double *T, const int i, const int j)
{
	double temp = T[i];
	T[i] = T[j];
	T[j] = temp;
}

inline void swap(double **T, const int i, const int j)
{
	double *temp = T[i];
	T[i] = T[j];
	T[j] = temp;
}

void centerMass(const int N, double **points, double *cm)
{
	for (size_t i = 0; i < 3 * N; i++) {
		cm[i] = 0;
	}
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < 3 * N; j++) {
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
			// swap(values, i, temp_right);
			// swap(points, temp_right, i);
			swap(values, i, temp_right);
			std::swap(points[temp_right], points[i]);
		}
	}

	swap(values, temp_right + 1, h);
	// swap(points, temp_right + 1, h);
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
	/*for (size_t i = 0; i < 3 * N; i++) {
		x_ref[i] = 0;
	}*/
	for (size_t i = 0; i < 3 * N; i++) {
		x_op[i] = (1 + r_op) * cm[i] - r_op * worst[i];
	}
}

void shrinkSimplex(const int N, double **points, double *values)
{
	for (size_t i = 1; i < N + 1; i++) {
		for (size_t j = 0; j < 3 * N; j++) {
			points[i][j] = points[0][j] - (points[i][j] - points[0][j]) * .5;
		}
	}
}

void freePoints(double **points, double *other, const int N)
{
	for (size_t i = 0; i < N + 1; i++) {
		if (points[i] == other) {
			// std::cout<<points[i]<<" "<<other<<"\n\n\n";
			// delete[] points[i];
			points[i] = nullptr;			
			// break;
		}
	}
}

void NelderMead(const int N, double **points)
{
	double *values = new double[N + 1]();
	double *cm = new double[3 * N]();
	double *x_inc = new double[3 * N]();
	double *x_ref = new double[3 * N]();
	double *x_exc = new double[3 * N]();
	double *x_exp = new double[3 * N]();
	
	for (size_t i = 0; i < N + 1; i++) {
		values[i] = ELJ(N, points[i]);
	}

	quicksort(points, values, 0, N);

	const double r_inc = -.5;
	const double r_exc = .5;
	const double r_ref = 1;
	const double r_exp = 2;

	int iterations = 0;
	const int maxIterations = 100000;
	const double real_val = -12.712062;
	const double acc = 0.001;
	double of_value = values[0];
	while (of_value > acc + real_val && iterations < maxIterations) {
		centerMass(N, points, cm);
		generate_point(x_ref, points[N], cm, N, r_ref);
		double x_ref_val = ELJ(N, x_ref);
		if (x_ref_val >= values[0] && x_ref_val < values[N - 1]) {
			delete[] points[N];
			points[N] = x_ref;
			values[N] = x_ref_val;
		}
		else if (x_ref_val < values[0]) {
			// Eigen::MatrixXd x_exp = (1 + r_exp) * (*cm) - r_exp * points[N];
			generate_point(x_exp, points[N], cm, N, r_exp);
			double x_exp_val = ELJ(N, x_exp);
			if (x_exp_val < x_ref_val) {
				delete[] points[N];
				points[N] = x_exp;
				values[N] = x_exp_val;
			}
			else {
				delete[] points[N];
				points[N] = x_ref;
				values[N] = x_ref_val;
			}
		}
		else if (x_ref_val >= values[N - 1] && x_ref_val < values[N]) {
			// Eigen::MatrixXd x_exc = (1 + r_exc) * (*cm) - r_exc * points[N];
			generate_point(x_exc, points[N], cm, N, r_exc);
			double x_exc_val = ELJ(N, x_exc);
			if (x_exc_val <= x_ref_val) {
				delete[] points[N];
				points[N] = x_exc;
				values[N] = x_exc_val;
			}
		}
		else if (values[N] <= x_ref_val) {
			// Eigen::MatrixXd x_inc = (1 + r_inc) * (*cm) - r_inc * points[N];
			generate_point(x_inc, points[N], cm, N, r_inc);
			double x_inc_val = ELJ(N, x_inc);
			if (x_inc_val < values[N]) {
				delete[] points[N];
				points[N] = x_inc;
				values[N] = x_inc_val;
			}
		}
		else {
			shrinkSimplex(N, points, values);
		}
		quicksort(points, values, 0, N);

		of_value = values[0];
		std::cout << iterations << ": " << of_value << "\n";
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
	std::default_random_engine generator;
	std::uniform_real_distribution<double> dist(-2.5, 2.5);

	// Eigen::MatrixXd cluster = Eigen::MatrixXd(N, 3);
	// cluster << 0.7430002202, 0.2647603899, -0.0468575389, -0.7430002647, -0.2647604843, 0.0468569750, 0.1977276118, -0.4447220146,
	// 0.6224700350, 	-0.1977281310, 0.4447221826, -0.6224697723, -0.1822009635, 0.5970484122, 0.4844363476, 0.1822015272, -0.5970484858,
	// -0.4844360463;

	double **p = new double *[N + 1]; //{0.7430002202, 0.2647603899, -0.0468575389, -0.7430002647, -0.2647604843, 0.0468569750, 0.1977276118,
									  //-0.4447220146, 0.6224700350,
									  // -0.1977281310, 0.4447221826, -0.6224697723, -0.1822009635, 0.5970484122, 0.4844363476, 0.1822015272,
									  // -0.5970484858, -0.4844360463};

	/*Eigen::MatrixXd *points = new Eigen::MatrixXd[N + 1];
	for (int i = 0; i < N + 1; i++) {
		points[i] = Eigen::MatrixXd::Random(N, 3);
	}
	for (size_t t = 0; t < N + 1; t++) {
		int count = 0;
		p[t] = new double[3 * N];
		for (size_t z = 0; z < N; z++) {
			for (size_t i = 0; i < 3; i++) {
				p[t][count++] = points[t](z, i);
			}
		}
	}*/
	/*auto t1 = high_resolution_clock::now();
	std::cout << "ELJ with Eigen " << ELJ(N, points[0]) << "\n";
	auto t2 = high_resolution_clock::now();
	std::cout << "Time " << duration_cast<milliseconds>(t2 - t1).count() << "\n";
	t1 = high_resolution_clock::now();
	std::cout << "ELJ with array " << ELJ(N, p) << "\n";
	t2 = high_resolution_clock::now();
	std::cout << "Time " << duration_cast<milliseconds>(t2 - t1).count() << "\n";*/
	for (size_t t = 0; t < N + 1; t++) {
		p[t] = new double[3 * N];
		for (size_t z = 0; z < 3 * N; z++) {
			p[t][z] = dist(generator);
		}
	}

	NelderMead(N, p);

	// delete[] points;
	for (size_t t = 0; t < N + 1; t++) {
		if (p[t] != nullptr) {
			delete[] p[t];
			p[t] = nullptr;
		}
	}
	delete[] p;
	return 0;
}

/*
cluster << 0.7417086996	,	-0.94012601	,	1.6794889913	,
1.1428995281	,	0.1309386561	,	-1.6267780141	,
0.5821815101	,	0.1602065612	,	1.8664613577	,
1.2224824756	,	-1.11332137	,	0.6732462604	,
0.8681736347	,	1.411420051	,	0.1067586936	,
0.2613786608	,	-0.4998347644	,	-1.9036400963	,
-1.2605099808	,	-0.004005849	,	1.3399284071	,
-1.5456839794	,	0.2704608194	,	-0.8521687614	,
-0.4898758952	,	-1.5361885943	,	-0.0819993743	,
-0.4133918431	,	0.6332326434	,	1.6862968909	,
0.1231525509	,	0.60380389	,	-1.7714203435	,
-1.2643639318	,	1.1493250802	,	-0.2205544823	,
-0.7841514553	,	0.994222581	,	-1.2378809669	,
1.3910892599	,	0.4692592622	,	0.5175802218	,
-1.1177919632	,	1.0106135829	,	0.8845714552	,
0.985931256	,	-0.8777326463	,	-1.1374560893	,
1.2405629447	,	0.6141504884	,	-0.631489021	,
0.6106543229	,	-1.3948923191	,	-0.2163705381	,
-0.2628130067	,	-0.4755134653	,	1.512962229	,
0.1239871635	,	-1.2430230601	,	0.7994458456	,
0.446684728	,	0.8745291261	,	1.0265148964	,
-0.6435922675	,	-0.1135401305	,	-1.3809850985	,
-1.405218175	,	0.1358438271	,	0.2451530457	,
-0.1133757753	,	-1.0203793834	,	-0.9934342167	,
-0.8757560583	,	-0.7702015206	,	0.6295555589	,
0.2215734281	,	1.0908117928	,	-0.7823141263	,
-1.0169880498	,	-0.6321319558	,	-0.4667580573	,
-0.2573941217	,	1.2463003872	,	0.2288539788	,
1.0842547623	,	-0.3914965837	,	-0.1533306775	,
0.4285236487	,	0.4391324441	,	0.0453073312	,
0.5915269789	,	-0.2377665171	,	0.8497816929	,
0.3438160493	,	-0.0233378728	,	-0.9070424105	,
-0.3983522718	,	0.2193946681	,	0.680522414	,
-0.541083938	,	0.3546019478	,	-0.3933536098	,
-0.0202388891	,	-0.5347557666	,	-0.0154533868	;*/