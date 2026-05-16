#ifndef FITTING
#define FITTING
/*!
*      \file Fitting.h
*      \brief EM Algorithm
*	   \author Wei Feng
*      \date 04/28/2026
*
*/
#include "Algo.h"
#include <chrono>
#include "Logger.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

namespace acpd
{
	namespace fit2d
	{


		inline CoeffSet CreateCoeffSet2D(int deg)
		{
			if (deg < 0) {
				throw std::invalid_argument("CreateCoeffSet2D: deg must be nonnegative.");
			}

			CoeffSet coeffs(deg + 1);

			for (int j = 0; j <= deg; ++j) {
				coeffs[j] = MatrixXd::Zero(2, j + 1);
			}

			return coeffs;
		}


		inline void GetMatB4Analytic(
			MatrixXd& matB,
			int deg,
			const Eigen::Ref<const VectorXd>& nomiVal,
			const Eigen::Ref<const MatrixXd>& point0)
		{
			const int point_num = static_cast<int>(point0.rows());
			const int sn = Sn(2, 2, deg + 1);

			if (point0.cols() < 2) {
				throw std::invalid_argument("GetMatB4AnalyticEigen: point0 must have at least 2 columns.");
			}
			if (nomiVal.size() < 2) {
				throw std::invalid_argument("GetMatB4AnalyticEigen: nomiVal must have at least 2 entries.");
			}

			matB.setZero(2 * point_num, sn);

			for (int i = 0; i < point_num; ++i)
			{
				int offset = 0;

				const double dx = point0(i, 0) - nomiVal(0);
				const double dy = point0(i, 1) - nomiVal(1);

				for (int j = 0; j <= deg; ++j)
				{
					const double one_fac = 1.0 / Factorial(j);

					for (int k = 0; k <= j; ++k)
					{
						const int comb_num = Combination(j, k);
						const double basis =
							one_fac *
							static_cast<double>(comb_num) *
							std::pow(dx, j - k) *
							std::pow(dy, k);

						// coefficient block for the x component
						matB(2 * i + 0, offset + k) = basis;

						// coefficient block for the y component
						matB(2 * i + 1, offset + (j + 1) + k) = basis;
					}

					offset += 2 * (j + 1);
				}
			}
		}

		inline void GetL4Analytic(
			VectorXd& l,
			const CoeffSet& mat_set,
			int deg,
			const Eigen::Ref<const VectorXd>& nomiVal,
			const Eigen::Ref<const MatrixXd>& point0,
			const Eigen::Ref<const MatrixXd>& point)
		{
			ValidateCoeffSet(mat_set, 2, deg, "GetL4AnalyticEigen");

			const int point_num = static_cast<int>(point0.rows());

			if (point.rows() != point_num) {
				throw std::invalid_argument("GetL4AnalyticEigen: point0 and point row mismatch.");
			}
			if (point0.cols() < 2 || point.cols() < 2) {
				throw std::invalid_argument("GetL4AnalyticEigen: point0 and point must have at least 2 columns.");
			}
			if (nomiVal.size() < 2) {
				throw std::invalid_argument("GetL4AnalyticEigen: nomiVal must have at least 2 entries.");
			}

			l.setZero(2 * point_num);

			for (int i = 0; i < point_num; ++i)
			{
				l(2 * i + 0) = point(i, 0);
				l(2 * i + 1) = point(i, 1);

				const double dx = point0(i, 0) - nomiVal(0);
				const double dy = point0(i, 1) - nomiVal(1);

				for (int j = 0; j <= deg; ++j)
				{
					const double one_fac = 1.0 / Factorial(j);

					for (int k = 0; k <= j; ++k)
					{
						const int comb_num = Combination(j, k);
						const double basis =
							one_fac *
							static_cast<double>(comb_num) *
							std::pow(dx, j - k) *
							std::pow(dy, k);

						// estimated x component of the current mapping
						l(2 * i + 0) -= basis * mat_set[j](0, k);

						// estimated y component of the current mapping
						l(2 * i + 1) -= basis * mat_set[j](1, k);
					}
				}
			}
		}


		inline double GetDeltaVec4Analytic(
			CoeffSet& delta_set,
			const CoeffSet& mat_set,
			const Eigen::Ref<const VectorXd>& nomiVal,
			int deg,
			const Eigen::Ref<const MatrixXd>& point0,
			const Eigen::Ref<const MatrixXd>& point,
			const Eigen::Ref<const MatrixXd>& W)
		{
			ValidateCoeffSet(mat_set, 2, deg, "GetDeltaVec4AnalyticEigen.mat_set");

			if (point0.rows() != point.rows()) {
				throw std::invalid_argument("GetDeltaVec4AnalyticEigen: point0 and point row mismatch.");
			}
			if (point0.cols() < 2 || point.cols() < 2) {
				throw std::invalid_argument("GetDeltaVec4AnalyticEigen: point0 and point must have at least 2 columns.");
			}

			const int point_num = static_cast<int>(point0.rows());
			const int sn = Sn(2, 2, deg + 1);

			if (W.rows() != 2 * point_num || W.cols() != 2 * point_num) {
				throw std::invalid_argument("GetDeltaVec4AnalyticEigen: W must be 2M x 2M.");
			}

			MatrixXd matB = MatrixXd::Zero(2 * point_num, sn);
			VectorXd l = VectorXd::Zero(2 * point_num);
			VectorXd deltaV = VectorXd::Zero(sn);

			GetMatB4Analytic(matB, deg, nomiVal, point0);
			GetL4Analytic(l, mat_set, deg, nomiVal, point0, point);

			const bool ok = GetCorrectionWeighted(deltaV, matB, l, W);
			if (!ok) {
				throw std::runtime_error("GetDeltaVec4AnalyticEigen: weighted correction failed.");
			}

			delta_set = CreateCoeffSet2D(deg);

			for (int j = 0; j <= deg; ++j)
			{
				const int offset = Sn(2, 2, j);

				for (int k = 0; k <= j; ++k)
				{
					delta_set[j](0, k) = deltaV(offset + k);
					delta_set[j](1, k) = deltaV(offset + (j + 1) + k);
				}
			}

			return deltaV.norm();
		}

		inline double AnalyticFitting(
			CoeffSet& mat_set,
			const Eigen::Ref<const VectorXd>& nomiVal,
			int deg,
			const Eigen::Ref<const MatrixXd>& point0,
			const Eigen::Ref<const MatrixXd>& point,
			const Eigen::Ref<const MatrixXd>& W,
			double eps,
			int jt)
		{
			ValidateCoeffSet(mat_set, 2, deg, "AnalyticFittingEigen.mat_set");

			if (point0.rows() != point.rows()) {
				throw std::invalid_argument("AnalyticFittingEigen: point0 and pointAff row mismatch.");
			}
			if (point0.cols() < 2 || point.cols() < 2) {
				throw std::invalid_argument("AnalyticFittingEigen: point0 and pointAff must have at least 2 columns.");
			}
			if (W.rows() != 2 * point0.rows() || W.cols() != 2 * point0.rows()) {
				throw std::invalid_argument("AnalyticFittingEigen: W must be 2M x 2M.");
			}

			double mat_error = 0.0;

			for (int iter = 0; iter < jt; ++iter)
			{
				CoeffSet delta_set;

				mat_error = GetDeltaVec4Analytic(
					delta_set,
					mat_set,
					nomiVal,
					deg,
					point0,
					point,
					W
					);

				for (int j = 0; j <= deg; ++j)
				{
					for (int k = 0; k <= j; ++k)
					{
						mat_set[j](0, k) += delta_set[j](0, k);
						mat_set[j](1, k) += delta_set[j](1, k);
					}
				}

				if (mat_error < eps)
				{
					break;
				}
			}

			return mat_error;
		}

		inline double CalcT(VectorXd &X0, VectorXd &hk, MatrixXd &matB, VectorXd &l)
		{
			VectorXd bn = l - matB*X0;
			VectorXd fh0 = matB*hk;

			double a2 = (fh0.transpose()*fh0)(0);
			double ab = (fh0.transpose()*bn)(0);
			return ab / a2;
		}

		inline MatrixXd ApplyAnalyticMap(
			const CoeffSet& mat_set,
			const Eigen::Ref<const VectorXd>& nomiVal,
			int deg,
			const Eigen::Ref<const MatrixXd>& points)
		{
			ValidateCoeffSet(mat_set, 2, deg, "ApplyAnalyticMapEigen.mat_set");

			const int point_num = static_cast<int>(points.rows());

			if (points.cols() < 2) {
				throw std::invalid_argument("ApplyAnalyticMapEigen: points must have at least 2 columns.");
			}
			if (nomiVal.size() < 2) {
				throw std::invalid_argument("ApplyAnalyticMapEigen: nomiVal must have at least 2 entries.");
			}

			MatrixXd out = MatrixXd::Zero(point_num, 2);

			for (int i = 0; i < point_num; ++i)
			{
				const double dx = points(i, 0) - nomiVal(0);
				const double dy = points(i, 1) - nomiVal(1);

				double vx = 0.0;
				double vy = 0.0;

				for (int j = 0; j <= deg; ++j)
				{
					const double one_fac = 1.0 / Factorial(j);

					for (int k = 0; k <= j; ++k)
					{
						const int comb_num = Combination(j, k);
						const double basis =
							one_fac *
							static_cast<double>(comb_num) *
							std::pow(dx, j - k) *
							std::pow(dy, k);

						vx += basis * mat_set[j](0, k);
						vy += basis * mat_set[j](1, k);
					}
				}

				out(i, 0) = vx;
				out(i, 1) = vy;
			}

			return out;
		}

		inline MatrixXd BuildWeightMatrix2D(
			const Eigen::Ref<const VectorXd>& pointWeights)
		{
			const int M = static_cast<int>(pointWeights.size());

			MatrixXd W = MatrixXd::Zero(2 * M, 2 * M);

			for (int i = 0; i < M; ++i)
			{
				const double wi = pointWeights(i);
				W(2 * i + 0, 2 * i + 0) = wi;
				W(2 * i + 1, 2 * i + 1) = wi;
			}

			return W;
		}

		inline void ResizeCoeffSet2D(
			CoeffSet& coeffs,
			int deg,
			bool preserveOld)
		{
			if (deg < 0) {
				throw std::invalid_argument("ResizeCoeffSet2D: deg must be nonnegative.");
			}

			CoeffSet newCoeffs(deg + 1);

			for (int j = 0; j <= deg; ++j) {
				newCoeffs[j] = Eigen::MatrixXd::Zero(2, j + 1);
			}

			if (preserveOld) {
				const int oldMax = static_cast<int>(coeffs.size()) - 1;
				const int copyMax = min(oldMax, deg);

				for (int j = 0; j <= copyMax; ++j) {
					if (coeffs[j].rows() >= 2 && coeffs[j].cols() >= j + 1) {
						newCoeffs[j].block(0, 0, 2, j + 1) =
							coeffs[j].block(0, 0, 2, j + 1);
					}
				}
			}

			coeffs.swap(newCoeffs);
		}
		inline double AMVFF2D_MStep(
			Eigen::MatrixXd& Ycur,                         // M x 2, input current point set and output updated point set
			const Eigen::Ref<const Eigen::MatrixXd>& Zfit, // M_eff x 2
			const Eigen::Ref<const Eigen::MatrixXd>& Yfit, // M_eff x 2
			const Eigen::Ref<const Eigen::VectorXd>& rho,  // M_eff
			const Eigen::Ref<const Eigen::VectorXd>& nomiVal,
			int deg,
			double fitting_eps = 1e-10,
			int innerIter = 1)
		{
			if (Yfit.rows() != Zfit.rows()) {
				throw std::invalid_argument("AMVFF2D_MStep: Yfit and Zfit row mismatch.");
			}
			if (rho.size() != Yfit.rows()) {
				throw std::invalid_argument("AMVFF2D_MStep: rho size mismatch.");
			}
			if (Ycur.cols() < 2 || Yfit.cols() < 2 || Zfit.cols() < 2) {
				throw std::invalid_argument("AMVFF2D_MStep: all point matrices must have at least 2 columns.");
			}

			CoeffSet mat_set = CreateCoeffSet2D(deg);

			// Directly fit Yfit -> Zfit, so this round can start from zero coefficients.
			// If an identity-like incremental mapping is desired, use InitializeIdentity2D(mat_set, deg) instead.
			SetCoeffSetZero(mat_set);

			Eigen::MatrixXd W = BuildWeightMatrix2D(rho);

			const double err = AnalyticFitting(
				mat_set,
				nomiVal,
				deg,
				Yfit,
				Zfit,
				W,
				fitting_eps,
				innerIter
				);

			Eigen::MatrixXd Ynew = ApplyAnalyticMap(
				mat_set,
				nomiVal,
				deg,
				Ycur
				);

			Ycur = std::move(Ynew);

			return err;
		}


		inline double Analytic_CPD_2D(
			Vector2d* point_set4regist,
			int num4regist,
			Vector2d* aim_point_set,
			int aim_num,
			int cpdIterCount,
			double eps,
			int maxDeg = 10,
			double w = 0.0,
			const std::string& logFile = "",
			bool detailedLog = false)
		{
			if (point_set4regist == nullptr || aim_point_set == nullptr) {
				throw std::invalid_argument("Analytic_CPD_2D: null point array.");
			}
			if (num4regist <= 0 || aim_num <= 0) {
				throw std::invalid_argument("Analytic_CPD_2D: invalid point number.");
			}
			if (cpdIterCount <= 0) {
				throw std::invalid_argument("Analytic_CPD_2D: cpdIterCount must be positive.");
			}
			if (maxDeg < 1) {
				throw std::invalid_argument("Analytic_CPD_2D: maxDeg should be at least 1.");
			}
			if (!(w >= 0.0 && w < 1.0)) {
				throw std::invalid_argument("Analytic_CPD_2D: w must satisfy 0 <= w < 1.");
			}

			constexpr int d = 2;

			const double eps_sigma = 1e-12;
			const double eps_rho = 1e-12;

			/*
			* External RMSE rebound guard.
			* It is meaningful only when fixed and moving point sets have valid
			* pointwise correspondence, which is usually true for synthetic tests.
			*/
			const int externalRisePatience = 1;
			const double externalRelRiseTol = 0.02;   // 2% above best observed RMSE
			const double externalAbsRiseTol = 1e-12;

			/*
			* Internal patience parameters.
			*/
			const int stablePatience = 5;
			const int noImprovePatience = 8;
			const int minIterBeforeStop = 5;

			const double relImproveTol = 1e-6;
			const double absImproveTol = 1e-12;
			const double relRiseTol = 1e-3;
			const double absRiseTol = 1e-12;

			auto tic = std::chrono::high_resolution_clock::now();

			RegistrationCSVLogger logger;
			if (!logFile.empty()) {
				logger.open(logFile, detailedLog, true);
			}

			// fixed target X: N x 2
			Eigen::MatrixXd X = Vector2dArrayToMatrix(aim_point_set, aim_num);

			// current moving Y: M x 2
			Eigen::MatrixXd Ycur = Vector2dArrayToMatrix(point_set4regist, num4regist);

			const bool hasPointwiseMetric = (X.rows() == Ycur.rows());

			// fixed expansion center c = 0
			Eigen::VectorXd nomiVal(2);
			nomiVal << 0.0, 0.0;

			// sigma^2 initialization
			double sigma2 = InitializeSigma2(X, Ycur, eps_sigma);

			const double initialSoftScore =
				std::sqrt(static_cast<double>(d) * sigma2);

			double initialExternalRMSE =
				std::numeric_limits<double>::quiet_NaN();

			if (hasPointwiseMetric) {
				initialExternalRMSE = GetRegistRMSE(X, Ycur);
			}

			/*
			* Internal best state.
			* This uses the CPD-style soft score and does not require pointwise
			* correspondence.
			*/
			Eigen::MatrixXd YbestInternal = Ycur;

			double bestInternalScore = initialSoftScore;
			double bestInternalSoft = initialSoftScore;
			double bestInternalRMSE = initialExternalRMSE;
			int bestInternalIter = -1;

			/*
			* External best state.
			* This is used only when pointwise RMSE is meaningful.
			*/
			Eigen::MatrixXd YbestExternal = Ycur;

			double bestExternalRMSE = initialExternalRMSE;
			double bestExternalSoft = initialSoftScore;
			int bestExternalIter = -1;

			int externalRiseCount = 0;

			double prevInternalScore = bestInternalScore;

			int noImproveCount = 0;
			int stableCount = 0;
			int prevDeg = -1;

			double finalScore = hasPointwiseMetric ? initialExternalRMSE : initialSoftScore;

			for (int iter = 0; iter < cpdIterCount; ++iter)
			{
				// ------------------------------------------------------------
				// 1. E-step: posterior correspondence estimation
				// ------------------------------------------------------------
				PosteriorResult post = ComputePosteriorP(
					X,
					Ycur,
					sigma2,
					w,
					eps_sigma
					);

				if (post.NP <= eps_rho) {
					break;
				}

				// ------------------------------------------------------------
				// 2. Posterior condensation: soft targets and weights
				// ------------------------------------------------------------
				WeightedFitData fitData = BuildWeightedFitDataFromPosterior(
					Ycur,
					post,
					eps_rho,
					true
					);

				if (fitData.Yfit.rows() < 3) {
					// 2D affine fitting needs enough effective points.
					break;
				}

				// ------------------------------------------------------------
				// 3. Degree schedule and feasibility check
				// ------------------------------------------------------------
				const int rawDeg =
					DegreeScheduleDecreasingStages(iter, cpdIterCount, maxDeg);

				int deg = rawDeg;

				while (deg > 1 &&
					2 * fitData.Yfit.rows() < TotalParamCount(2, deg))
				{
					--deg;
				}

				if (2 * fitData.Yfit.rows() < TotalParamCount(2, deg)) {
					break;
				}

				const bool degreeChanged = (prevDeg >= 0 && deg != prevDeg);

				if (degreeChanged) {
					noImproveCount = 0;
					stableCount = 0;
					externalRiseCount = 0;
				}

				prevDeg = deg;

				Eigen::MatrixXd Yold = Ycur;

				// ------------------------------------------------------------
				// 4. M-step: weighted 2D structured analytic mapping fitting
				// ------------------------------------------------------------
				const double fitErr = AMVFF2D_MStep(
					Ycur,
					fitData.Zfit,
					fitData.Yfit,
					fitData.weights,
					nomiVal,
					deg,
					1e-10,
					1
					);

				// ------------------------------------------------------------
				// 5. Update sigma^2
				// ------------------------------------------------------------
				const double sigma2Next = UpdateSigma2FromStats(
					X,
					Ycur,
					post.rowSums,
					post.colSums,
					post.PX,
					eps_sigma
					);

				const double softScore =
					std::sqrt(static_cast<double>(d) * sigma2Next);

				double externalRMSE =
					std::numeric_limits<double>::quiet_NaN();

				if (hasPointwiseMetric) {
					externalRMSE = GetRegistRMSE(X, Ycur);
				}

				finalScore = hasPointwiseMetric ? externalRMSE : softScore;

				const double deltaY =
					(Ycur - Yold).norm() / (Yold.norm() + 1e-12);

				const double deltaSigma =
					std::abs(sigma2Next - sigma2) /
					(std::abs(sigma2) + 1e-12);

				double deltaInternalScore =
					std::numeric_limits<double>::infinity();

				if (std::isfinite(prevInternalScore) && std::isfinite(softScore)) {
					deltaInternalScore =
						std::abs(prevInternalScore - softScore) /
						(std::abs(prevInternalScore) + 1e-12);
				}

				sigma2 = sigma2Next;

				// ------------------------------------------------------------
				// 6. Update internal best state
				// ------------------------------------------------------------
				bool internalImproved = false;

				if (std::isfinite(softScore)) {
					const double requiredImprove =
						max(absImproveTol,
							relImproveTol * std::abs(bestInternalScore));

					if (!std::isfinite(bestInternalScore) ||
						bestInternalScore - softScore > requiredImprove)
					{
						internalImproved = true;
					}
				}

				if (internalImproved) {
					bestInternalScore = softScore;
					bestInternalSoft = softScore;
					bestInternalRMSE = externalRMSE;
					YbestInternal = Ycur;
					bestInternalIter = iter;
					noImproveCount = 0;
				}
				else {
					++noImproveCount;
				}

				// ------------------------------------------------------------
				// 7. Update external best state and rebound counter
				// ------------------------------------------------------------
				if (hasPointwiseMetric && std::isfinite(externalRMSE)) {
					bool externalImproved = false;

					if (!std::isfinite(bestExternalRMSE)) {
						externalImproved = true;
					}
					else {
						const double requiredExternalImprove =
							max(absImproveTol,
								relImproveTol * std::abs(bestExternalRMSE));

						if (bestExternalRMSE - externalRMSE > requiredExternalImprove) {
							externalImproved = true;
						}
					}

					if (externalImproved) {
						bestExternalRMSE = externalRMSE;
						bestExternalSoft = softScore;
						YbestExternal = Ycur;
						bestExternalIter = iter;
						externalRiseCount = 0;
					}
					else {
						const bool externalClearlyWorse =
							std::isfinite(bestExternalRMSE) &&
							externalRMSE >
							bestExternalRMSE * (1.0 + externalRelRiseTol)
							+ externalAbsRiseTol;

						if (externalClearlyWorse) {
							++externalRiseCount;
						}
						else {
							externalRiseCount = 0;
						}
					}
				}

				// ------------------------------------------------------------
				// 8. Stable convergence counter
				// ------------------------------------------------------------
				if (deltaY < eps &&
					deltaSigma < eps &&
					deltaInternalScore < eps)
				{
					++stableCount;
				}
				else {
					stableCount = 0;
				}

				prevInternalScore = softScore;

				// ------------------------------------------------------------
				// 9. Logging
				// ------------------------------------------------------------
				auto toc_in = std::chrono::high_resolution_clock::now();

				const double runtimeSec =
					std::chrono::duration_cast<std::chrono::microseconds>(
						toc_in - tic
						).count() / 1000000.0;

				std::printf(
					"iter=%d, deg=%d, soft=%.15f, extRMSE=%.15f, "
					"bestExtRMSE=%.15f, extRise=%d, runtime=%.15f\n",
					iter,
					deg,
					softScore,
					externalRMSE,
					bestExternalRMSE,
					externalRiseCount,
					runtimeSec
					);

				if (logger.is_open()) {
					RegistrationIterLog rec;
					rec.iter = iter;
					rec.degree = deg;
					rec.rmse = externalRMSE;
					rec.sigma2 = sigma2Next;
					rec.fitErr = fitErr;
					rec.runtimeSec = runtimeSec;
					rec.activeCount = static_cast<int>(fitData.Yfit.rows());
					rec.NP = post.NP;

					logger.write(rec);
				}

				// ------------------------------------------------------------
				// 10. Stopping tests
				// ------------------------------------------------------------

				// Internal soft residual is already tiny.
				if (softScore < eps) {
					break;
				}

				/*
				* External RMSE rebound guard.
				* Mainly for synthetic experiments with known pointwise correspondence.
				*/
				if (hasPointwiseMetric &&
					iter >= minIterBeforeStop &&
					externalRiseCount >= externalRisePatience)
				{
					break;
				}

				/*
				* Internal convergence tests.
				*/
				if (iter >= minIterBeforeStop) {
					if (stableCount >= stablePatience) {
						break;
					}

					const bool internalClearlyWorse =
						std::isfinite(softScore) &&
						std::isfinite(bestInternalScore) &&
						softScore >
						bestInternalScore * (1.0 + relRiseTol) + absRiseTol;

					if (noImproveCount >= stablePatience && internalClearlyWorse) {
						break;
					}

					if (noImproveCount >= noImprovePatience) {
						break;
					}
				}
			}

			// ------------------------------------------------------------
			// 11. Roll back to the best state
			// ------------------------------------------------------------
			if (hasPointwiseMetric && std::isfinite(bestExternalRMSE)) {
				Ycur = YbestExternal;
				finalScore = bestExternalRMSE;
			}
			else if (std::isfinite(bestInternalScore)) {
				Ycur = YbestInternal;
				finalScore = bestInternalScore;
			}

			// copy result back
			MatrixToVector2dArray(Ycur, point_set4regist, num4regist);

			return finalScore;
		}
	}

	namespace fit3d
	{
		// ============================================================
		// Weighted fit data for 3D M-step
		// ============================================================

		struct WeightedFitData3D
		{
			MatrixXd Yfit;       // M_eff x 3
			MatrixXd Zfit;       // M_eff x 3
			VectorXd weights;    // M_eff
			std::vector<int> activeIndex;
		};

		inline WeightedFitData3D BuildWeightedFitData3DFromPosterior(
			const Eigen::Ref<const MatrixXd>& Ycur, // M x 3
			const PosteriorResult& post,
			double eps_rho = 1e-12,
			bool dropTinyWeightRows = true)
		{
			const int M = static_cast<int>(Ycur.rows());

			if (Ycur.cols() < 3) {
				throw std::invalid_argument("BuildWeightedFitData3DFromPosterior: Ycur must have at least 3 columns.");
			}
			if (post.rowSums.size() != M) {
				throw std::invalid_argument("BuildWeightedFitData3DFromPosterior: rowSums size mismatch.");
			}
			if (post.PX.rows() != M || post.PX.cols() < 3) {
				throw std::invalid_argument("BuildWeightedFitData3DFromPosterior: PX size mismatch.");
			}

			WeightedFitData3D out;
			out.activeIndex.reserve(M);

			for (int m = 0; m < M; ++m)
			{
				if (!dropTinyWeightRows || post.rowSums(m) > eps_rho) {
					out.activeIndex.push_back(m);
				}
			}

			const int Meff = static_cast<int>(out.activeIndex.size());

			out.Yfit = MatrixXd::Zero(Meff, 3);
			out.Zfit = MatrixXd::Zero(Meff, 3);
			out.weights = VectorXd::Zero(Meff);

			for (int i = 0; i < Meff; ++i)
			{
				const int m = out.activeIndex[i];
				const double rho = max(post.rowSums(m), eps_rho);

				out.Yfit.row(i) = Ycur.row(m).leftCols(3);
				out.Zfit.row(i) = post.PX.row(m).leftCols(3) / rho;
				out.weights(i) = rho;
			}

			return out;
		}

		// ============================================================
		// Array conversion
		// ============================================================

		inline MatrixXd Vector3dArrayToMatrix(
			const Eigen::Vector3d* pts,
			int num)
		{
			if (pts == nullptr || num <= 0) {
				throw std::invalid_argument("Vector3dArrayToMatrix: invalid input.");
			}

			MatrixXd M(num, 3);

			for (int i = 0; i < num; ++i)
			{
				M(i, 0) = pts[i](0);
				M(i, 1) = pts[i](1);
				M(i, 2) = pts[i](2);
			}

			return M;
		}

		inline void MatrixToVector3dArray(
			const Eigen::Ref<const MatrixXd>& M,
			Eigen::Vector3d* pts,
			int num)
		{
			if (pts == nullptr || num <= 0) {
				throw std::invalid_argument("MatrixToVector3dArray: invalid output.");
			}
			if (M.rows() != num || M.cols() < 3) {
				throw std::invalid_argument("MatrixToVector3dArray: matrix size mismatch.");
			}

			for (int i = 0; i < num; ++i)
			{
				pts[i](0) = M(i, 0);
				pts[i](1) = M(i, 1);
				pts[i](2) = M(i, 2);
			}
		}

		// Optional external RMSE if point order is known.
		inline double GetRegistRMSE3D(
			const Eigen::Ref<const MatrixXd>& fixed,
			const Eigen::Ref<const MatrixXd>& re)
		{
			if (fixed.rows() != re.rows() || fixed.cols() < 3 || re.cols() < 3) {
				throw std::invalid_argument("GetRegistRMSE3D: size mismatch.");
			}

			if (fixed.rows() <= 0) {
				return 0.0;
			}

			const MatrixXd diff = fixed.leftCols(3) - re.leftCols(3);

			return std::sqrt(diff.squaredNorm() / static_cast<double>(fixed.rows()));
		}

		// ============================================================
		// One M-step: fit Yfit -> Zfit and update all Ycur
		// ============================================================

		inline double AMVFF3D_MStep(
			MatrixXd& Ycur,                                 // M x 3, updated in-place
			const Eigen::Ref<const MatrixXd>& Yfit,          // M_eff x 3
			const Eigen::Ref<const MatrixXd>& Zfit,          // M_eff x 3
			const Eigen::Ref<const VectorXd>& weights,       // M_eff
			const Eigen::Ref<const VectorXd>& nomiVal,       // 3
			int deg,
			double fitting_eps = 1e-10,
			int innerIter = 1)
		{
			if (Yfit.rows() != Zfit.rows()) {
				throw std::invalid_argument("AMVFF3D_MStep: Yfit and Zfit row mismatch.");
			}
			if (weights.size() != Yfit.rows()) {
				throw std::invalid_argument("AMVFF3D_MStep: weights size mismatch.");
			}

			CoeffSet mat_set = CreateCoeffSet(3, deg);

			// Direct fit Yfit -> Zfit.
			SetCoeffSetZero(mat_set);

			const double err = AnalyticFitting3D(
				mat_set,
				nomiVal,
				deg,
				Yfit,
				Zfit,
				weights,
				fitting_eps,
				innerIter
				);

			MatrixXd Ynew = ApplyAnalyticMap3D(
				mat_set,
				nomiVal,
				deg,
				Ycur
				);

			Ycur = std::move(Ynew);

			return err;
		}

		// ============================================================
		// Analytic-CPD 3D main loop
		// ============================================================

		inline double Analytic_CPD_3D(
			Eigen::Vector3d* point_set4regist,
			int num4regist,
			Eigen::Vector3d* aim_point_set,
			int aim_num,
			int cpdIterCount,
			double eps,
			int maxDeg = 10,
			double w = 0.1,
			const std::string& logFile = "",
			bool detailedLog = false)
		{
			if (point_set4regist == nullptr || aim_point_set == nullptr) {
				throw std::invalid_argument("Analytic_CPD_3D: null point array.");
			}
			if (num4regist <= 0 || aim_num <= 0) {
				throw std::invalid_argument("Analytic_CPD_3D: invalid point number.");
			}
			if (cpdIterCount <= 0) {
				throw std::invalid_argument("Analytic_CPD_3D: cpdIterCount must be positive.");
			}
			if (maxDeg < 1) {
				throw std::invalid_argument("Analytic_CPD_3D: maxDeg must be at least 1.");
			}
			if (!(w >= 0.0 && w < 1.0)) {
				throw std::invalid_argument("Analytic_CPD_3D: w must satisfy 0 <= w < 1.");
			}

			constexpr int d = 3;

			const double eps_sigma = 1e-12;
			const double eps_rho = 1e-12;

			/*
			* This rebound guard uses the external pointwise RMSE.
			* It is meaningful only for synthetic or ordered datasets where X and Y
			* have one-to-one pointwise correspondence.
			*
			* For real unordered point sets, this guard is automatically disabled
			* because X.rows() != Y.rows() or pointwise RMSE is not meaningful.
			*/
			const int externalRisePatience = 1;
			const double externalRelRiseTol = 0.02;   // 2% higher than the best observed RMSE
			const double externalAbsRiseTol = 1e-12;

			/*
			* Internal patience parameters.
			* These are used for soft-score / sigma / motion stabilization.
			*/
			const int stablePatience = 5;
			const int noImprovePatience = 8;
			const int minIterBeforeStop = 5;

			const double relImproveTol = 1e-6;
			const double absImproveTol = 1e-12;
			const double relRiseTol = 1e-3;
			const double absRiseTol = 1e-12;

			auto tic = std::chrono::high_resolution_clock::now();

			RegistrationCSVLogger logger;
			if (!logFile.empty()) {
				logger.open(logFile, detailedLog, true);
			}

			// fixed target X: N x 3
			Eigen::MatrixXd X = Vector3dArrayToMatrix(aim_point_set, aim_num);

			// current moving Y: M x 3
			Eigen::MatrixXd Ycur = Vector3dArrayToMatrix(point_set4regist, num4regist);

			const bool hasPointwiseMetric = (X.rows() == Ycur.rows());

			Eigen::VectorXd nomiVal = Eigen::VectorXd::Zero(3);

			double sigma2 = InitializeSigma2(X, Ycur, eps_sigma);

			const double initialSoftScore =
				std::sqrt(static_cast<double>(d) * sigma2);

			double initialExternalRMSE =
				std::numeric_limits<double>::quiet_NaN();

			if (hasPointwiseMetric) {
				initialExternalRMSE = GetRegistRMSE3D(X, Ycur);
			}

			/*
			* Internal best state.
			* This uses the CPD-style soft score and does not require pointwise
			* correspondence.
			*/
			Eigen::MatrixXd YbestInternal = Ycur;

			double bestInternalScore = initialSoftScore;
			double bestInternalSoft = initialSoftScore;
			double bestInternalRMSE = initialExternalRMSE;
			int bestInternalIter = -1;

			/*
			* External best state.
			* This is used only when pointwise RMSE is meaningful.
			* It is especially useful in synthetic experiments to prevent returning
			* a later overfitted state after RMSE rebound.
			*/
			Eigen::MatrixXd YbestExternal = Ycur;

			double bestExternalRMSE = initialExternalRMSE;
			double bestExternalSoft = initialSoftScore;
			int bestExternalIter = -1;

			int externalRiseCount = 0;

			double prevInternalScore = bestInternalScore;

			int noImproveCount = 0;
			int stableCount = 0;
			int prevDeg = -1;

			double finalScore = hasPointwiseMetric ? initialExternalRMSE : initialSoftScore;

			for (int iter = 0; iter < cpdIterCount; ++iter)
			{
				// ------------------------------------------------------------
				// 1. E-step: posterior correspondence estimation
				// ------------------------------------------------------------
				PosteriorResult post = ComputePosteriorP(
					X,
					Ycur,
					sigma2,
					w,
					eps_sigma
					);

				if (post.NP <= eps_rho) {
					break;
				}

				// ------------------------------------------------------------
				// 2. Posterior condensation: soft targets and weights
				// ------------------------------------------------------------
				WeightedFitData3D fitData = BuildWeightedFitData3DFromPosterior(
					Ycur,
					post,
					eps_rho,
					true
					);

				if (fitData.Yfit.rows() < 4) {
					// 3D affine fitting needs enough effective points.
					break;
				}

				// ------------------------------------------------------------
				// 3. Degree schedule and feasibility check
				// ------------------------------------------------------------
				const int rawDeg =
					DegreeScheduleDecreasingStages(iter, cpdIterCount, maxDeg);

				int deg = rawDeg;

				while (deg > 1 &&
					3 * fitData.Yfit.rows() < TotalParamCount(3, deg))
				{
					--deg;
				}

				if (3 * fitData.Yfit.rows() < TotalParamCount(3, deg)) {
					break;
				}

				const bool degreeChanged = (prevDeg >= 0 && deg != prevDeg);

				if (degreeChanged) {
					noImproveCount = 0;
					stableCount = 0;
					externalRiseCount = 0;
				}

				prevDeg = deg;

				Eigen::MatrixXd Yold = Ycur;

				// ------------------------------------------------------------
				// 4. M-step: weighted 3D structured analytic mapping fitting
				// ------------------------------------------------------------
				const double fitErr = AMVFF3D_MStep(
					Ycur,
					fitData.Yfit,
					fitData.Zfit,
					fitData.weights,
					nomiVal,
					deg,
					1e-10,
					1
					);

				// ------------------------------------------------------------
				// 5. Update sigma^2
				// ------------------------------------------------------------
				const double sigma2Next = UpdateSigma2FromStats(
					X,
					Ycur,
					post.rowSums,
					post.colSums,
					post.PX,
					eps_sigma
					);

				const double softScore =
					std::sqrt(static_cast<double>(d) * sigma2Next);

				double externalRMSE =
					std::numeric_limits<double>::quiet_NaN();

				if (hasPointwiseMetric) {
					externalRMSE = GetRegistRMSE3D(X, Ycur);
				}

				finalScore = hasPointwiseMetric ? externalRMSE : softScore;

				const double deltaY =
					(Ycur - Yold).norm() / (Yold.norm() + 1e-12);

				const double deltaSigma =
					std::abs(sigma2Next - sigma2) /
					(std::abs(sigma2) + 1e-12);

				double deltaInternalScore =
					std::numeric_limits<double>::infinity();

				if (std::isfinite(prevInternalScore) && std::isfinite(softScore)) {
					deltaInternalScore =
						std::abs(prevInternalScore - softScore) /
						(std::abs(prevInternalScore) + 1e-12);
				}

				sigma2 = sigma2Next;

				// ------------------------------------------------------------
				// 6. Update internal best state
				// ------------------------------------------------------------
				bool internalImproved = false;

				if (std::isfinite(softScore)) {
					const double requiredImprove =
						max(absImproveTol,
							relImproveTol * std::abs(bestInternalScore));

					if (!std::isfinite(bestInternalScore) ||
						bestInternalScore - softScore > requiredImprove)
					{
						internalImproved = true;
					}
				}

				if (internalImproved) {
					bestInternalScore = softScore;
					bestInternalSoft = softScore;
					bestInternalRMSE = externalRMSE;
					YbestInternal = Ycur;
					bestInternalIter = iter;
					noImproveCount = 0;
				}
				else {
					++noImproveCount;
				}

				// ------------------------------------------------------------
				// 7. Update external best state and rebound counter
				// ------------------------------------------------------------
				if (hasPointwiseMetric && std::isfinite(externalRMSE)) {
					bool externalImproved = false;

					if (!std::isfinite(bestExternalRMSE)) {
						externalImproved = true;
					}
					else {
						const double requiredExternalImprove =
							max(absImproveTol,
								relImproveTol * std::abs(bestExternalRMSE));

						if (bestExternalRMSE - externalRMSE > requiredExternalImprove) {
							externalImproved = true;
						}
					}

					if (externalImproved) {
						bestExternalRMSE = externalRMSE;
						bestExternalSoft = softScore;
						YbestExternal = Ycur;
						bestExternalIter = iter;
						externalRiseCount = 0;
					}
					else {
						const bool externalClearlyWorse =
							std::isfinite(bestExternalRMSE) &&
							externalRMSE >
							bestExternalRMSE * (1.0 + externalRelRiseTol)
							+ externalAbsRiseTol;

						if (externalClearlyWorse) {
							++externalRiseCount;
						}
						else {
							externalRiseCount = 0;
						}
					}
				}

				// ------------------------------------------------------------
				// 8. Stable convergence counter
				// ------------------------------------------------------------
				if (deltaY < eps &&
					deltaSigma < eps &&
					deltaInternalScore < eps)
				{
					++stableCount;
				}
				else {
					stableCount = 0;
				}

				prevInternalScore = softScore;

				// ------------------------------------------------------------
				// 9. Logging
				// ------------------------------------------------------------
				auto toc_in = std::chrono::high_resolution_clock::now();

				const double runtimeSec =
					std::chrono::duration_cast<std::chrono::microseconds>(
						toc_in - tic
						).count() / 1000000.0;

				std::printf(
					"iter=%d, deg=%d, soft=%.15f, extRMSE=%.15f, "
					"bestExtRMSE=%.15f, extRise=%d, runtime=%.15f\n",
					iter,
					deg,
					softScore,
					externalRMSE,
					bestExternalRMSE,
					externalRiseCount,
					runtimeSec
					);

				if (logger.is_open()) {
					RegistrationIterLog rec;
					rec.iter = iter;
					rec.degree = deg;
					rec.rmse = externalRMSE;
					rec.sigma2 = sigma2Next;
					rec.fitErr = fitErr;
					rec.runtimeSec = runtimeSec;
					rec.activeCount = static_cast<int>(fitData.Yfit.rows());
					rec.NP = post.NP;

					logger.write(rec);
				}

				// ------------------------------------------------------------
				// 10. Stopping tests
				// ------------------------------------------------------------

				// Internal soft residual is already tiny.
				if (softScore < eps) {
					break;
				}

				/*
				* External RMSE rebound guard.
				* This is mainly for synthetic experiments with known pointwise
				* correspondence. If RMSE has moved clearly above the historical best
				* for several consecutive iterations, we stop and later roll back to
				* the best external state.
				*/
				if (hasPointwiseMetric &&
					iter >= minIterBeforeStop &&
					externalRiseCount >= externalRisePatience)
				{
					break;
				}

				/*
				* Internal convergence tests.
				* These are less dependent on external ground-truth ordering.
				*/
				if (iter >= minIterBeforeStop) {
					if (stableCount >= stablePatience) {
						break;
					}

					const bool internalClearlyWorse =
						std::isfinite(softScore) &&
						std::isfinite(bestInternalScore) &&
						softScore >
						bestInternalScore * (1.0 + relRiseTol) + absRiseTol;

					if (noImproveCount >= stablePatience && internalClearlyWorse) {
						break;
					}

					if (noImproveCount >= noImprovePatience) {
						break;
					}
				}
			}

			// ------------------------------------------------------------
			// 11. Roll back to the best state
			// ------------------------------------------------------------
			if (hasPointwiseMetric && std::isfinite(bestExternalRMSE)) {
				Ycur = YbestExternal;
				finalScore = bestExternalRMSE;
			}
			else if (std::isfinite(bestInternalScore)) {
				Ycur = YbestInternal;
				finalScore = bestInternalScore;
			}

			MatrixToVector3dArray(Ycur, point_set4regist, num4regist);

			return finalScore;
		}
	}
}


#endif