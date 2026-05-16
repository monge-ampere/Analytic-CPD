#pragma once
/*!
*      \file Algo.h
*      \brief algorithm
*	   \author Wei Feng
*      \date 04/28/2026
*/
#ifndef MATHEMATICAL
#define MATHEMATICAL
#include<math.h>
#include <stdlib.h>
#include <time.h> 
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <random>
#include <stdexcept>

using namespace Eigen;
using namespace std;

#define INFINITYNUM 16777216
#define MAX_VERTEX_NUM 32768

#define PI 3.14159265358979323846264338327950288419716939937510


using CoeffSet = std::vector<Eigen::MatrixXd>;

struct PosteriorResult
{
	Eigen::MatrixXd P;        // M x N
	Eigen::VectorXd rowSums;  // M x 1, rho = P * 1
	Eigen::VectorXd colSums;  // N x 1, P^T * 1
	Eigen::MatrixXd PX;       // M x d
	double NP = 0.0;          // sum_{m,n} p_mn
	double c = 0.0;           // uniform outlier term
};

struct WeightedFitData
{
	Eigen::MatrixXd Yfit;       // M_eff x 2
	Eigen::MatrixXd Zfit;       // M_eff x 2
	Eigen::VectorXd weights;    // M_eff

	vector<int> activeIndex;
};


#define RADTODEG(ang)			((ang) * 180.0 / PI)

/*!
* the point is moved by a randomly generated analytic map.
* \param v target
* \param mo the moved point
* \param num the point number
* \param deg the order of the Taylor series
*/
template<class TT1>
void __declspec(dllexport) AnalyticTransf2d(TT1(*v)[2], TT1(*mo)[2], int num, int deg)
{
	srand((unsigned)time(NULL));
	int sn = Sn(2, 2, deg + 1);
	double *pCoef = new double[sn];
	double nomiVal[2] = { 0 };

	FILE *fp = fopen("E:\\coef.csv", "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != sn; ++i)
	{
		pCoef[i] = 0.8 * rand() / RAND_MAX - 0.4;

		fprintf(fp, "%f\n", pCoef[i]);

	}

	nomiVal[0] = 0.1;
	nomiVal[1] = 0.1;

	fclose(fp);

	pCoef[2] += 1;
	pCoef[5] += 1;

	double *pTerm = new double[2];
	int sum = 0;

	for (int i = 0; i != num; ++i)
	{
		mo[i][0] = 0;
		mo[i][1] = 0;
		for (int j = 1; j <= deg; ++j)
		{
			sum = Sn(2, 2, j);
			for (int k = 0; k <= j; ++k)
			{
				pTerm[0] = pCoef[sum + k * 2 + 0]
					* Combination(j, k)
					*pow(v[i][0] - nomiVal[0], j - k)
					*pow(v[i][1] - nomiVal[1], k) / Factorial(j);

				pTerm[1] = pCoef[sum + k * 2 + 1]
					* Combination(j, k)
					*pow(v[i][0] - nomiVal[0], j - k)
					*pow(v[i][1] - nomiVal[1], k) / Factorial(j);

				mo[i][0] += pTerm[0];
				mo[i][1] += pTerm[1];
			}
		}

	}

	delete[]pTerm;
	delete[]pCoef;
}

inline bool GetCorrectionWeighted(
	VectorXd& correction,
	const Eigen::Ref<const MatrixXd>& matrixB,
	const Eigen::Ref<const VectorXd>& l,
	const Eigen::Ref<const MatrixXd>& W)
{
	if (matrixB.rows() != l.size()) {
		throw std::invalid_argument("GetCorrectionWeightedEigen: matrixB and l size mismatch.");
	}
	if (W.rows() != matrixB.rows() || W.cols() != matrixB.rows()) {
		throw std::invalid_argument("GetCorrectionWeightedEigen: W size mismatch.");
	}

	const MatrixXd nbb = matrixB.transpose() * W * matrixB;
	const VectorXd rhs = matrixB.transpose() * W * l;

	Eigen::LDLT<MatrixXd> solver(nbb);
	if (solver.info() != Eigen::Success) {
		return false;
	}

	correction = solver.solve(rhs);

	if (solver.info() != Eigen::Success || !correction.allFinite()) {
		return false;
	}

	return true;
}

inline int DegreeScheduleDecreasingStages(
	int iter,
	int maxIter,
	int maxDeg)
{
	if (maxDeg < 1) {
		throw std::invalid_argument(
			"DegreeScheduleDecreasingStages: maxDeg must be >= 1.");
	}

	if (maxIter <= 0) {
		return maxDeg;
	}

	if (iter >= maxIter) {
		return maxDeg;
	}

	if (iter < 0) {
		iter = 0;
	}

	const int D = maxDeg;

	const long long totalWeight =
		static_cast<long long>(D) * static_cast<long long>(D + 1) / 2;

	const long long fullUnits = maxIter / totalWeight;
	const int remainder = static_cast<int>(maxIter % totalWeight);

	long long cumulative = 0;

	for (int degree = 1; degree <= D; ++degree) {
		long long duration =
			fullUnits * static_cast<long long>(D - degree + 1);

		int rest = remainder;

		for (int prefix = D; prefix >= 1 && rest > 0; --prefix) {
			const int addCount = min(prefix, rest);

			if (degree <= addCount) {
				++duration;
			}

			rest -= addCount;
		}

		cumulative += duration;

		if (iter < cumulative) {
			return degree;
		}
	}

	return D;
}


/*!
* the sum of n terms of the arithmetic sequence
*/
template<class T>
T __declspec(dllexport) Sn(T a0, T d, int n)
{
	T v = n*a0 + n*(n - 1)*d / 2;
	return v;
}

/*!
* the sum of n terms of the 1,3,6,10^ sequence
*/
inline int Sn136(int n)
{
	int v = n*(n + 1)*(n + 2) / 2;
	return v;
}

inline double Factorial(int n)
{
	if (n < 0) {
		throw std::invalid_argument("Factorial: n must be nonnegative.");
	}

	double v = 1.0;
	for (int i = 2; i <= n; ++i) {
		v *= static_cast<double>(i);
	}
	return v;
}

inline int Combination(int n, int k)
{
	if (k < 0 || k > n) {
		return 0;
	}
	if (k == 0 || k == n) {
		return 1;
	}

	k = min(k, n - k);

	long long result = 1;
	for (int i = 1; i <= k; ++i) {
		result = result * (n - k + i) / i;
	}

	return static_cast<int>(result);
}

// Update sigma^2 given rowSums = P.rowwise().sum(),
// colSums = P.colwise().sum().transpose(), and PX = P * X.
inline double UpdateSigma2FromStats(
	const Eigen::Ref<const Eigen::MatrixXd>& X,        // N x d
	const Eigen::Ref<const Eigen::MatrixXd>& Ynew,     // M x d
	const Eigen::Ref<const Eigen::VectorXd>& rowSums,  // M
	const Eigen::Ref<const Eigen::VectorXd>& colSums,  // N
	const Eigen::Ref<const Eigen::MatrixXd>& PX,       // M x d
	double eps_sigma = 1e-12)
{
	const int N = static_cast<int>(X.rows());
	const int d = static_cast<int>(X.cols());
	const int M = static_cast<int>(Ynew.rows());

	if (Ynew.cols() != d) {
		throw std::invalid_argument("UpdateSigma2FromStats: X and Ynew must have the same point dimension.");
	}
	if (rowSums.size() != M || colSums.size() != N) {
		throw std::invalid_argument("UpdateSigma2FromStats: rowSums/colSums dimension mismatch.");
	}
	if (PX.rows() != M || PX.cols() != d) {
		throw std::invalid_argument("UpdateSigma2FromStats: PX must be of size M x d.");
	}

	const double NP = rowSums.sum();
	if (NP <= eps_sigma) {
		return eps_sigma;
	}

	const double term1 = colSums.dot(X.rowwise().squaredNorm());
	const double term2 = (PX.array() * Ynew.array()).sum();
	const double term3 = rowSums.dot(Ynew.rowwise().squaredNorm());

	double sigma2 = (term1 - 2.0 * term2 + term3) / (NP * static_cast<double>(d));

	if (!std::isfinite(sigma2) || sigma2 < eps_sigma) {
		sigma2 = eps_sigma;
	}

	return sigma2;
}

inline double InitializeSigma2(
	const Eigen::Ref<const Eigen::MatrixXd>& X,   // N x d
	const Eigen::Ref<const Eigen::MatrixXd>& Y,   // M x d
	double eps_sigma = 1e-12)
{
	const int N = static_cast<int>(X.rows());
	const int d = static_cast<int>(X.cols());
	const int M = static_cast<int>(Y.rows());

	if (Y.cols() != d) {
		throw std::invalid_argument("InitializeSigma2: X and Y must have the same point dimension.");
	}
	if (d <= 0 || M <= 0 || N <= 0) {
		throw std::invalid_argument("InitializeSigma2: invalid empty input.");
	}

	double acc = 0.0;
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < M; ++m) {
			acc += (X.row(n) - Y.row(m)).squaredNorm();
		}
	}

	double sigma2 = acc / (static_cast<double>(d) * M * N);
	if (!std::isfinite(sigma2) || sigma2 < eps_sigma) {
		sigma2 = eps_sigma;
	}
	return sigma2;
}

inline PosteriorResult ComputePosteriorP(
	const Eigen::Ref<const Eigen::MatrixXd>& X,   // N x d
	const Eigen::Ref<const Eigen::MatrixXd>& Y,   // M x d
	double sigma2,
	double w,
	double eps_sigma = 1e-12)
{
	const int N = static_cast<int>(X.rows());
	const int d = static_cast<int>(X.cols());
	const int M = static_cast<int>(Y.rows());

	if (Y.cols() != d) {
		throw std::invalid_argument("ComputePosteriorP: X and Y must have the same point dimension.");
	}
	if (N <= 0 || M <= 0 || d <= 0) {
		throw std::invalid_argument("ComputePosteriorP: empty input.");
	}
	if (!(w >= 0.0 && w < 1.0)) {
		throw std::invalid_argument("ComputePosteriorP: w must satisfy 0 <= w < 1.");
	}

	sigma2 = max(sigma2, eps_sigma);

	PosteriorResult out;
	out.P = Eigen::MatrixXd::Zero(M, N);
	out.rowSums = Eigen::VectorXd::Zero(M);
	out.colSums = Eigen::VectorXd::Zero(N);
	out.PX = Eigen::MatrixXd::Zero(M, d);

	// c = (2*pi*sigma2)^(d/2) * w/(1-w) * M/N
	if (w == 0.0) {
		out.c = 0.0;
	}
	else {
		out.c = std::pow(2.0 * PI * sigma2, 0.5 * static_cast<double>(d))
			* (w / (1.0 - w))
			* (static_cast<double>(M) / static_cast<double>(N));
	}

	const double invTwoSigma2 = 1.0 / (2.0 * sigma2);

	// Column-wise computation: for each fixed x_n, compute the posterior over all y_m.
	Eigen::VectorXd kcol(M);

	for (int n = 0; n < N; ++n) {
		double denom = out.c;

		// First compute all Gaussian affinities in this column.
		for (int m = 0; m < M; ++m) {
			const double dist2 = (X.row(n) - Y.row(m)).squaredNorm();
			const double kval = std::exp(-dist2 * invTwoSigma2);
			kcol[m] = kval;
			denom += kval;
		}

		if (denom <= eps_sigma) {
			denom = eps_sigma;
		}

		// Normalize to obtain the posterior, while accumulating the required statistics.
		for (int m = 0; m < M; ++m) {
			const double p = kcol[m] / denom;
			out.P(m, n) = p;
			out.rowSums[m] += p;
			out.colSums[n] += p;
			out.PX.row(m).noalias() += p * X.row(n);
		}
	}

	out.NP = out.rowSums.sum();
	return out;
}

inline WeightedFitData BuildWeightedFitDataFromPosterior(
	const Eigen::Ref<const Eigen::MatrixXd>& Ycur,   // M x 2
	const PosteriorResult& post,
	double eps_rho = 1e-12,
	bool dropTinyWeightRows = true)
{
	const int M = static_cast<int>(Ycur.rows());
	const int d = static_cast<int>(Ycur.cols());

	if (d < 2) {
		throw std::invalid_argument("BuildWeightedFitDataFromPosterior: Ycur must have at least 2 columns.");
	}
	if (post.rowSums.size() != M) {
		throw std::invalid_argument("BuildWeightedFitDataFromPosterior: rowSums size mismatch.");
	}
	if (post.PX.rows() != M || post.PX.cols() != d) {
		throw std::invalid_argument("BuildWeightedFitDataFromPosterior: PX size mismatch.");
	}

	WeightedFitData out;
	out.activeIndex.reserve(M);

	for (int m = 0; m < M; ++m) {
		if (!dropTinyWeightRows || post.rowSums[m] > eps_rho) {
			out.activeIndex.push_back(m);
		}
	}

	const int Meff = static_cast<int>(out.activeIndex.size());

	out.Yfit = Eigen::MatrixXd::Zero(Meff, 2);
	out.Zfit = Eigen::MatrixXd::Zero(Meff, 2);
	out.weights = Eigen::VectorXd::Zero(Meff);

	for (int i = 0; i < Meff; ++i) {
		const int m = out.activeIndex[i];
		const double rho = max(post.rowSums[m], eps_rho);

		out.Yfit.row(i) = Ycur.row(m).leftCols(2);
		out.Zfit.row(i) = post.PX.row(m).leftCols(2) / rho;
		out.weights(i) = rho;
	}

	return out;
}

inline Eigen::MatrixXd Vector2dArrayToMatrix(
	const Eigen::Vector2d* pts,
	int num)
{
	if (pts == nullptr || num <= 0) {
		throw std::invalid_argument("Vector2dArrayToMatrix: invalid input.");
	}

	Eigen::MatrixXd M(num, 2);

	for (int i = 0; i < num; ++i) {
		M(i, 0) = pts[i](0);
		M(i, 1) = pts[i](1);
	}

	return M;
}

inline void MatrixToVector2dArray(
	const Eigen::Ref<const Eigen::MatrixXd>& M,
	Eigen::Vector2d* pts,
	int num)
{
	if (pts == nullptr || num <= 0) {
		throw std::invalid_argument("MatrixToVector2dArray: invalid output.");
	}
	if (M.rows() != num || M.cols() < 2) {
		throw std::invalid_argument("MatrixToVector2dArray: matrix size mismatch.");
	}

	for (int i = 0; i < num; ++i) {
		pts[i](0) = M(i, 0);
		pts[i](1) = M(i, 1);
	}
}

inline double GetRegistRMSE(
	const Eigen::Ref<const Eigen::MatrixXd>& fixed,
	const Eigen::Ref<const Eigen::MatrixXd>& re)
{
	if (fixed.rows() != re.rows() || fixed.cols() != re.cols()) {
		throw std::invalid_argument("GetRegistRMSE: size mismatch.");
	}

	const int num = static_cast<int>(fixed.rows());

	if (num <= 0) {
		return 0.0;
	}

	const Eigen::MatrixXd diff = fixed - re;

	return std::sqrt(diff.squaredNorm() / static_cast<double>(num));
}


// ============================================================
// Basic helpers
// ============================================================

inline int MonomialCount(int dim, int order)
{
	if (dim <= 0 || order < 0) {
		throw std::invalid_argument("MonomialCount: invalid dim or order.");
	}

	return Combination(order + dim - 1, dim - 1);
}

inline int TotalMonomialCount(int dim, int deg)
{
	if (dim <= 0 || deg < 0) {
		throw std::invalid_argument("TotalMonomialCount: invalid dim or deg.");
	}

	return Combination(deg + dim, dim);
}

inline int TotalParamCount(int dim, int deg)
{
	return dim * TotalMonomialCount(dim, deg);
}

// ============================================================
// Coefficient-set utilities
// ============================================================

inline CoeffSet CreateCoeffSet(int dim, int deg)
{
	if (dim <= 0 || deg < 0) {
		throw std::invalid_argument("CreateCoeffSet: invalid dim or deg.");
	}

	CoeffSet coeffs(deg + 1);

	for (int j = 0; j <= deg; ++j) {
		const int cnum = MonomialCount(dim, j);
		coeffs[j] = Eigen::MatrixXd::Zero(dim, cnum);
	}

	return coeffs;
}

inline void ValidateCoeffSet(
	const CoeffSet& coeffs,
	int dim,
	int deg,
	const char* name)
{
	if (dim <= 0 || deg < 0) {
		throw std::invalid_argument(std::string(name) + ": invalid dim or deg.");
	}

	if (static_cast<int>(coeffs.size()) < deg + 1) {
		throw std::invalid_argument(std::string(name) + ": coeffs size is too small.");
	}

	for (int j = 0; j <= deg; ++j) {
		const int cnum = MonomialCount(dim, j);

		if (coeffs[j].rows() < dim || coeffs[j].cols() < cnum) {
			throw std::invalid_argument(std::string(name) + ": coefficient block has wrong size.");
		}
	}
}

inline void SetCoeffSetZero(CoeffSet& coeffs)
{
	for (auto& C : coeffs) {
		C.setZero();
	}
}

inline void InitializeIdentity(CoeffSet& mat_set, int dim, int deg)
{
	ValidateCoeffSet(mat_set, dim, deg, "InitializeIdentity");

	SetCoeffSetZero(mat_set);

	if (deg >= 1) {
		for (int i = 0; i < dim; ++i) {
			mat_set[1](i, i) = 1.0;
		}
	}
}


// ============================================================
// Taylor coefficient generation
// ============================================================
//
// This follows your original GetCoefOfTaylor order:
//
// for a = order down to 0:
//     rem = order - a
//     for c = 0 to rem:
//         b = rem - c
//
// monomial = x^a y^b z^c / (a! b! c!)
//
// ============================================================

inline VectorXd GetTaylorCoef3D(
	const Eigen::Ref<const VectorXd>& p,
	const Eigen::Ref<const VectorXd>& nomiVal,
	int order)
{
	if (p.size() < 3 || nomiVal.size() < 3) {
		throw std::invalid_argument("GetTaylorCoef3D: p and nomiVal must have at least 3 entries.");
	}

	const double x = p(0) - nomiVal(0);
	const double y = p(1) - nomiVal(1);
	const double z = p(2) - nomiVal(2);

	const int cnum = MonomialCount(3, order);
	VectorXd coef = VectorXd::Zero(cnum);

	int index = 0;

	for (int a = order; a >= 0; --a)
	{
		const int rem = order - a;

		for (int c = 0; c <= rem; ++c)
		{
			const int b = rem - c;

			const double v =
				std::pow(x, a) *
				std::pow(y, b) *
				std::pow(z, c) /
				(Factorial(a) * Factorial(b) * Factorial(c));

			coef(index++) = v;
		}
	}

	return coef;
}

// ============================================================
// Build design matrix B
// ============================================================

inline void GetMatB4Analytic3D(
	MatrixXd& matB,
	int deg,
	const Eigen::Ref<const VectorXd>& nomiVal,
	const Eigen::Ref<const MatrixXd>& point0) // M x 3
{
	if (point0.cols() < 3) {
		throw std::invalid_argument("GetMatB4Analytic3DEigen: point0 must have at least 3 columns.");
	}
	if (nomiVal.size() < 3) {
		throw std::invalid_argument("GetMatB4Analytic3DEigen: nomiVal must have at least 3 entries.");
	}

	const int point_num = static_cast<int>(point0.rows());
	const int param_num = TotalParamCount(3, deg);

	matB.setZero(3 * point_num, param_num);

	for (int i = 0; i < point_num; ++i)
	{
		VectorXd p(3);
		p << point0(i, 0), point0(i, 1), point0(i, 2);

		int offset = 0;

		for (int j = 0; j <= deg; ++j)
		{
			const VectorXd coef = GetTaylorCoef3D(p, nomiVal, j);
			const int cnum = static_cast<int>(coef.size());

			for (int k = 0; k < cnum; ++k)
			{
				// x output block
				matB(3 * i + 0, offset + k) = coef(k);

				// y output block
				matB(3 * i + 1, offset + cnum + k) = coef(k);

				// z output block
				matB(3 * i + 2, offset + 2 * cnum + k) = coef(k);
			}

			offset += 3 * cnum;
		}
	}
}

// ============================================================
// Build residual vector l
// ============================================================

inline void GetL4Analytic3D(
	VectorXd& l,
	const CoeffSet& mat_set,
	int deg,
	const Eigen::Ref<const VectorXd>& nomiVal,
	const Eigen::Ref<const MatrixXd>& point0, // M x 3
	const Eigen::Ref<const MatrixXd>& point)  // M x 3 target
{
	ValidateCoeffSet(mat_set, 3, deg, "GetL4Analytic3DEigen.mat_set");

	if (point0.rows() != point.rows()) {
		throw std::invalid_argument("GetL4Analytic3DEigen: point0 and point row mismatch.");
	}
	if (point0.cols() < 3 || point.cols() < 3) {
		throw std::invalid_argument("GetL4Analytic3DEigen: point0 and point must have at least 3 columns.");
	}
	if (nomiVal.size() < 3) {
		throw std::invalid_argument("GetL4Analytic3DEigen: nomiVal must have at least 3 entries.");
	}

	const int point_num = static_cast<int>(point0.rows());

	l.setZero(3 * point_num);

	for (int i = 0; i < point_num; ++i)
	{
		l(3 * i + 0) = point(i, 0);
		l(3 * i + 1) = point(i, 1);
		l(3 * i + 2) = point(i, 2);

		VectorXd p(3);
		p << point0(i, 0), point0(i, 1), point0(i, 2);

		for (int j = 0; j <= deg; ++j)
		{
			const VectorXd coef = GetTaylorCoef3D(p, nomiVal, j);
			const int cnum = static_cast<int>(coef.size());

			for (int k = 0; k < cnum; ++k)
			{
				l(3 * i + 0) -= coef(k) * mat_set[j](0, k);
				l(3 * i + 1) -= coef(k) * mat_set[j](1, k);
				l(3 * i + 2) -= coef(k) * mat_set[j](2, k);
			}
		}
	}
}

// ============================================================
// Weighted least-squares correction
// Recommended: use point weights directly instead of explicit W.
// Each point weight rho_i is applied to its x,y,z equations.
// ============================================================

inline bool GetCorrectionWeightedByPoint3D(
	VectorXd& correction,
	const Eigen::Ref<const MatrixXd>& matrixB,
	const Eigen::Ref<const VectorXd>& l,
	const Eigen::Ref<const VectorXd>& pointWeights)
{
	if (matrixB.rows() != l.size()) {
		throw std::invalid_argument("GetCorrectionWeightedByPoint3D: matrixB and l size mismatch.");
	}

	const int obs = static_cast<int>(matrixB.rows());
	if (obs % 3 != 0) {
		throw std::invalid_argument("GetCorrectionWeightedByPoint3D: observation rows must be 3M.");
	}

	const int M = obs / 3;

	if (pointWeights.size() != M) {
		throw std::invalid_argument("GetCorrectionWeightedByPoint3D: pointWeights size mismatch.");
	}

	MatrixXd Bw = matrixB;
	VectorXd lw = l;

	for (int i = 0; i < M; ++i)
	{
		const double wi = max(pointWeights(i), 0.0);
		const double s = std::sqrt(wi);

		Bw.row(3 * i + 0) *= s;
		Bw.row(3 * i + 1) *= s;
		Bw.row(3 * i + 2) *= s;

		lw(3 * i + 0) *= s;
		lw(3 * i + 1) *= s;
		lw(3 * i + 2) *= s;
	}

	// More stable than solving normal equations explicitly.
	Eigen::ColPivHouseholderQR<MatrixXd> solver(Bw);

	correction = solver.solve(lw);

	return correction.allFinite();
}

// ============================================================
// One correction step
// ============================================================

inline double GetDeltaVec3D4Analytic(
	CoeffSet& delta_set,
	const CoeffSet& mat_set,
	const Eigen::Ref<const VectorXd>& nomiVal,
	int deg,
	const Eigen::Ref<const MatrixXd>& point0,       // source, M x 3
	const Eigen::Ref<const MatrixXd>& point,        // target, M x 3
	const Eigen::Ref<const VectorXd>& pointWeights) // M
{
	ValidateCoeffSet(mat_set, 3, deg, "GetDeltaVec3D4AnalyticEigen.mat_set");

	if (point0.rows() != point.rows()) {
		throw std::invalid_argument("GetDeltaVec3D4AnalyticEigen: point0 and point row mismatch.");
	}
	if (point0.cols() < 3 || point.cols() < 3) {
		throw std::invalid_argument("GetDeltaVec3D4AnalyticEigen: point0 and point must have at least 3 columns.");
	}
	if (pointWeights.size() != point0.rows()) {
		throw std::invalid_argument("GetDeltaVec3D4AnalyticEigen: pointWeights size mismatch.");
	}

	const int point_num = static_cast<int>(point0.rows());
	const int param_num = TotalParamCount(3, deg);

	if (3 * point_num < param_num) {
		// Not enough equations.
		throw std::runtime_error("GetDeltaVec3D4AnalyticEigen: underdetermined fitting system.");
	}

	MatrixXd matB = MatrixXd::Zero(3 * point_num, param_num);
	VectorXd l = VectorXd::Zero(3 * point_num);
	VectorXd deltaV = VectorXd::Zero(param_num);

	GetMatB4Analytic3D(matB, deg, nomiVal, point0);
	GetL4Analytic3D(l, mat_set, deg, nomiVal, point0, point);

	if (!matB.allFinite() || !l.allFinite()) {
		throw std::runtime_error("GetDeltaVec3D4AnalyticEigen: matB or l contains NaN/Inf.");
	}

	const bool ok = GetCorrectionWeightedByPoint3D(
		deltaV,
		matB,
		l,
		pointWeights
		);

	if (!ok) {
		throw std::runtime_error("GetDeltaVec3D4AnalyticEigen: weighted correction failed.");
	}

	delta_set = CreateCoeffSet(3, deg);

	int offset = 0;

	for (int j = 0; j <= deg; ++j)
	{
		const int cnum = MonomialCount(3, j);

		for (int k = 0; k < cnum; ++k)
		{
			delta_set[j](0, k) = deltaV(offset + k);
			delta_set[j](1, k) = deltaV(offset + cnum + k);
			delta_set[j](2, k) = deltaV(offset + 2 * cnum + k);
		}

		offset += 3 * cnum;
	}

	return deltaV.norm();
}

// ============================================================
// Analytic fitting, correction model
// ============================================================

inline double AnalyticFitting3D(
	CoeffSet& mat_set,
	const Eigen::Ref<const VectorXd>& nomiVal,
	int deg,
	const Eigen::Ref<const MatrixXd>& point0,       // source Yfit
	const Eigen::Ref<const MatrixXd>& pointAff,     // target Zfit
	const Eigen::Ref<const VectorXd>& pointWeights, // rho
	double eps = 1e-10,
	int jt = 1)
{
	ValidateCoeffSet(mat_set, 3, deg, "AnalyticFitting3DEigen.mat_set");

	if (point0.rows() != pointAff.rows()) {
		throw std::invalid_argument("AnalyticFitting3DEigen: point0 and pointAff row mismatch.");
	}
	if (pointWeights.size() != point0.rows()) {
		throw std::invalid_argument("AnalyticFitting3DEigen: pointWeights size mismatch.");
	}

	double mat_error = 0.0;

	for (int iter = 0; iter < jt; ++iter)
	{
		CoeffSet delta_set;

		mat_error = GetDeltaVec3D4Analytic(
			delta_set,
			mat_set,
			nomiVal,
			deg,
			point0,
			pointAff,
			pointWeights
			);

		for (int j = 0; j <= deg; ++j)
		{
			const int cnum = MonomialCount(3, j);

			for (int k = 0; k < cnum; ++k)
			{
				mat_set[j](0, k) += delta_set[j](0, k);
				mat_set[j](1, k) += delta_set[j](1, k);
				mat_set[j](2, k) += delta_set[j](2, k);
			}
		}

		if (mat_error < eps) {
			break;
		}
	}

	return mat_error;
}

// ============================================================
// Apply analytic mapping to point set
// ============================================================

inline MatrixXd ApplyAnalyticMap3D(
	const CoeffSet& mat_set,
	const Eigen::Ref<const VectorXd>& nomiVal,
	int deg,
	const Eigen::Ref<const MatrixXd>& points) // M x 3
{
	ValidateCoeffSet(mat_set, 3, deg, "ApplyAnalyticMap3DEigen.mat_set");

	if (points.cols() < 3) {
		throw std::invalid_argument("ApplyAnalyticMap3DEigen: points must have at least 3 columns.");
	}
	if (nomiVal.size() < 3) {
		throw std::invalid_argument("ApplyAnalyticMap3DEigen: nomiVal must have at least 3 entries.");
	}

	const int point_num = static_cast<int>(points.rows());

	MatrixXd out = MatrixXd::Zero(point_num, 3);

	for (int i = 0; i < point_num; ++i)
	{
		VectorXd p(3);
		p << points(i, 0), points(i, 1), points(i, 2);

		double vx = 0.0;
		double vy = 0.0;
		double vz = 0.0;

		for (int j = 0; j <= deg; ++j)
		{
			const VectorXd coef = GetTaylorCoef3D(p, nomiVal, j);
			const int cnum = static_cast<int>(coef.size());

			for (int k = 0; k < cnum; ++k)
			{
				vx += coef(k) * mat_set[j](0, k);
				vy += coef(k) * mat_set[j](1, k);
				vz += coef(k) * mat_set[j](2, k);
			}
		}

		out(i, 0) = vx;
		out(i, 1) = vy;
		out(i, 2) = vz;
	}

	return out;
}

inline Eigen::Vector3d AnalyticTrans3D(
	const CoeffSet& mat_set,
	int deg,
	const Eigen::Ref<const Eigen::Vector3d>& srcObj,
	const Eigen::Ref<const Eigen::Vector3d>& nomiVal)
{
	ValidateCoeffSet(mat_set, 3, deg, "AnalyticTrans3DEigen.mat_set");

	Eigen::Vector3d qstObj = Eigen::Vector3d::Zero();

	for (int order = 0; order <= deg; ++order)
	{
		Eigen::VectorXd coef = GetTaylorCoef3D(srcObj, nomiVal, order);
		const int cnum = static_cast<int>(coef.size());

		for (int k = 0; k < cnum; ++k)
		{
			qstObj(0) += mat_set[order](0, k) * coef(k);
			qstObj(1) += mat_set[order](1, k) * coef(k);
			qstObj(2) += mat_set[order](2, k) * coef(k);
		}
	}

	return qstObj;
}


inline double CompactBump(double r, double sigma)
{
	sigma = max(sigma, 1e-12);
	const double s = r / sigma;

	if (s >= 1.0) {
		return 0.0;
	}

	const double denom = 1.0 - s * s;
	return std::exp(-1.0 / denom);
}

inline void BumpBicentricWeights3D(
	const Eigen::Ref<const Eigen::Vector3d>& y,
	const Eigen::Ref<const Eigen::Vector3d>& c1,
	const Eigen::Ref<const Eigen::Vector3d>& c2,
	double sigma,
	double& w1,
	double& w2)
{
	const double b1 = CompactBump((y - c1).norm(), sigma);
	const double b2 = CompactBump((y - c2).norm(), sigma);

	const double s = b1 + b2;

	if (s <= 1e-16) {
		w1 = 0.5;
		w2 = 0.5;
	}
	else {
		w1 = b1 / s;
		w2 = b2 / s;
	}
}

inline void SaveCoeffSetToCSV3D(
	const CoeffSet& mat_set,
	int deg,
	const std::string& filename)
{
	ValidateCoeffSet(mat_set, 3, deg, "SaveCoeffSetToCSV3D.mat_set");

	std::ofstream fp(filename);
	if (!fp.is_open()) {
		throw std::runtime_error("SaveCoeffSetToCSV3D: failed to open file.");
	}

	for (int order = 0; order <= deg; ++order)
	{
		const int cnum = MonomialCount(3, order);

		for (int row = 0; row < 3; ++row)
		{
			for (int k = 0; k < cnum; ++k)
			{
				fp << order << "," << row << "," << k << "," << mat_set[order](row, k) << "\n";
			}
		}
	}
}

inline void AnalyticT3d(
	const Eigen::Ref<const Eigen::MatrixXd>& v, // num x 3, source / target prototype
	Eigen::MatrixXd& mo,                        // num x 3, moved output
	int deg,
	unsigned int seed = std::random_device{}(),
	bool saveCoef = false,
	const std::string& coefFile = "E:\\coef.csv",
	bool bicentricSecondOrder = true,
	bool cleanLinearIdentity = false)
{
	if (deg < 0) {
		throw std::invalid_argument("AnalyticT3dEigen: deg must be nonnegative.");
	}
	if (v.cols() < 3) {
		throw std::invalid_argument("AnalyticT3dEigen: v must have at least 3 columns.");
	}

	const int num = static_cast<int>(v.rows());
	if (num <= 0) {
		throw std::invalid_argument("AnalyticT3dEigen: empty point set.");
	}

	mo = Eigen::MatrixXd::Zero(num, 3);

	// fixed Taylor center, consistent with your old 3D generator
	Eigen::Vector3d nomiVal = Eigen::Vector3d::Zero();

	std::mt19937 rng(seed);

	std::uniform_real_distribution<double> baseDist(-0.4, 0.4);
	std::uniform_real_distribution<double> secondPerturbDist(-0.2, 0.2);

	// ------------------------------------------------------------
	// 1. Randomly initialize coefficients from order 0 to deg
	// ------------------------------------------------------------
	CoeffSet mat_set = CreateCoeffSet(3, deg);

	for (int order = 0; order <= deg; ++order)
	{
		const int cnum = MonomialCount(3, order);

		for (int r = 0; r < 3; ++r)
		{
			for (int k = 0; k < cnum; ++k)
			{
				mat_set[order](r, k) = baseDist(rng);
			}
		}
	}

	// ------------------------------------------------------------
	// 2. Set the first-order diagonal terms to 1
	// ------------------------------------------------------------
	if (deg >= 1)
	{
		if (cleanLinearIdentity)
		{
			mat_set[1].setZero();
		}

		// With the current 3D monomial order, order-1 basis is:
		// k=0: x
		// k=1: y
		// k=2: z
		mat_set[1](0, 0) = 1.0;
		mat_set[1](1, 1) = 1.0;
		mat_set[1](2, 2) = 1.0;
	}

	if (saveCoef)
	{
		SaveCoeffSetToCSV3D(mat_set, deg, coefFile);
	}

	// ------------------------------------------------------------
	// 3. If deg < 2 or bicentricSecondOrder disabled,
	//    use a single global analytic map directly.
	// ------------------------------------------------------------
	if (deg < 2 || !bicentricSecondOrder)
	{
		for (int i = 0; i < num; ++i)
		{
			Eigen::Vector3d p = v.row(i).leftCols(3).transpose();
			Eigen::Vector3d q = AnalyticTrans3D(mat_set, deg, p, nomiVal);

			mo(i, 0) = q(0);
			mo(i, 1) = q(1);
			mo(i, 2) = q(2);
		}

		return;
	}

	// ------------------------------------------------------------
	// 4. Bicentric second-order mapping
	// ------------------------------------------------------------

	// Choose two centers: min x and max x points
	int idx_minx = 0;
	int idx_maxx = 0;

	for (int i = 1; i < num; ++i)
	{
		if (v(i, 0) < v(idx_minx, 0)) {
			idx_minx = i;
		}
		if (v(i, 0) > v(idx_maxx, 0)) {
			idx_maxx = i;
		}
	}

	Eigen::Vector3d c1 = v.row(idx_minx).leftCols(3).transpose();
	Eigen::Vector3d c2 = v.row(idx_maxx).leftCols(3).transpose();

	const double d2c = (c1 - c2).squaredNorm();
	const double sigma2 = (d2c > 1e-12) ? (0.25 * d2c) : 1.0;
	const double sigma = std::sqrt(sigma2);

	// Backup original second-order coefficients
	Eigen::MatrixXd originalOrder2 = mat_set[2];

	// Generate two second-order coefficient blocks
	Eigen::MatrixXd Q1 = originalOrder2;
	Eigen::MatrixXd Q2 = originalOrder2;

	for (int r = 0; r < Q1.rows(); ++r)
	{
		for (int k = 0; k < Q1.cols(); ++k)
		{
			Q1(r, k) += secondPerturbDist(rng);
			Q2(r, k) += secondPerturbDist(rng);
		}
	}

	// ------------------------------------------------------------
	// 5. Pointwise evaluation:
	//    mat_set[2](y) = w1(y) Q1 + w2(y) Q2
	// ------------------------------------------------------------
	for (int i = 0; i < num; ++i)
	{
		Eigen::Vector3d p = v.row(i).leftCols(3).transpose();

		double w1 = 0.5;
		double w2 = 0.5;

		BumpBicentricWeights3D(p, c1, c2, sigma, w1, w2);

		mat_set[2] = w1 * Q1 + w2 * Q2;

		Eigen::Vector3d q = AnalyticTrans3D(mat_set, deg, p, nomiVal);

		mo(i, 0) = q(0);
		mo(i, 1) = q(1);
		mo(i, 2) = q(2);
	}

	// Restore optional; local variable anyway, so not necessary.
	mat_set[2] = originalOrder2;
}

inline void AnalyticT3dArrayWrapper(
	double(*v)[3],
	double(*mo)[3],
	int num,
	int deg,
	unsigned int seed)
{
	if (v == nullptr || mo == nullptr || num <= 0) {
		throw std::invalid_argument("AnalyticT3dArrayWrapper: invalid input.");
	}

	Eigen::MatrixXd Vin(num, 3);
	Eigen::MatrixXd Mout;

	for (int i = 0; i < num; ++i)
	{
		Vin(i, 0) = v[i][0];
		Vin(i, 1) = v[i][1];
		Vin(i, 2) = v[i][2];
	}

	AnalyticT3d(
		Vin,
		Mout,
		deg,
		seed,
		false
		);

	for (int i = 0; i < num; ++i)
	{
		mo[i][0] = Mout(i, 0);
		mo[i][1] = Mout(i, 1);
		mo[i][2] = Mout(i, 2);
	}
}
#endif