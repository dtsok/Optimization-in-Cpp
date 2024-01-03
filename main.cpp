#include <eigen3/Eigen/Dense>
#include <iostream>
#include <memory>

double ELJ(int N, const Eigen::MatrixXd &X, double epsilon = 1, double sigma = 1)
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

inline void swap(double *&T, const int i, const int j)
{
	double temp = T[i];
	T[i] = T[j];
	T[j] = temp;
}

void initialization(int N, Eigen::MatrixXd *&points)
{
	double min = ELJ(N, points[0]), max_1 = min, max_2 = min;
	int min_index = 0, max_1_index = 0, max_2_index = 0;
	double *values = new double[N + 1];
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

	delete[] values;
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

void NelderMead(int N, Eigen::MatrixXd *&points)
{
	initialization(N, points);

	std::unique_ptr<Eigen::MatrixXd> cm = centerMass(N, points);
}

int main(int argc, char const *argv[])
{
	int N = 1;
	if (argc > 1) {
		N = std::stoi(argv[1]);
	}

	// Eigen::MatrixXd cluster = Eigen::MatrixXd(N, 3);
	// cluster << 0.7430002202, 0.2647603899, -0.0468575389, -0.7430002647, -0.2647604843, 0.0468569750, 0.1977276118, -0.4447220146,
	// 0.6224700350, -0.1977281310, 0.4447221826, -0.6224697723, -0.1822009635, 0.5970484122, 0.4844363476, 0.1822015272, -0.5970484858,
	// -0.4844360463;
	Eigen::MatrixXd *points = new Eigen::MatrixXd[N + 1];
	for (int i = 0; i < N + 1; i++) {
		points[i] = Eigen::MatrixXd::Random(N, 3);
	}

	NelderMead(N, points);

	delete[] points;
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