// dllmain.cpp : 定义 DLL 应用程序的入口点。
#include "stdafx.h"
#include "Fitting.h"
#include <iostream>

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

inline double feature_nonrigid_regist2d(double(*mo)[2], int mnum
	, double(*fix)[2], int fnum, int iterCount)
{
	using namespace acpd;
	using namespace fit2d;
	Vector2d *moMat = new Vector2d[mnum];
	Vector2d *fixMat = new Vector2d[fnum];

	for (int i = 0; i != mnum; ++i)
	{
		moMat[i](0) = mo[i][0];
		moMat[i](1) = mo[i][1];
	}

	for (int i = 0; i != fnum; ++i)
	{
		fixMat[i](0) = fix[i][0];
		fixMat[i](1) = fix[i][1];
	}
	char strPara[MAX_PATH] = { 0 };
	int len = GetPrivateProfileStringA("Param",
		"Deg",
		"No Deg",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	int deg = atoi(strPara);

	len = GetPrivateProfileStringA("Param",
		"maxIterCount",
		"No maxIterCount",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	int maxIterCount = atoi(strPara);

	len = GetPrivateProfileStringA("Param",
		"Omega",
		"No Omega",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	double omega = atof(strPara);

	len = GetPrivateProfileStringA("Param",
		"Eps",
		"No Eps",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	double eps = atof(strPara);


	char error_time[MAX_PATH] = { 0 };

	len = GetPrivateProfileStringA("Param",
		"OutputErrorTime",
		"No OutputErrorTime",
		error_time,
		MAX_PATH,
		".//ACPD.ini");

	double error = 0;
	error = Analytic_CPD_2D(
		moMat,
		mnum,
		fixMat,
		fnum,
		maxIterCount,
		eps,
		deg,
		omega,
		error_time);

	for (int i = 0; i != mnum; ++i)
	{
		mo[i][0] = moMat[i](0);
		mo[i][1] = moMat[i](1);
	}

	delete[] moMat;
	delete[] fixMat;

	return error;
}

inline double feature_nonrigid_regist3d(double(*mo)[3], int mnum
	, double(*fix)[3], int fnum, int iterCount)
{
	using namespace acpd;
	using namespace fit3d;
	Vector3d *moMat = new Vector3d[mnum];
	Vector3d *fixMat = new Vector3d[fnum];

	for (int i = 0; i != mnum; ++i)
	{
		moMat[i](0) = mo[i][0];
		moMat[i](1) = mo[i][1];
		moMat[i](2) = mo[i][2];
	}

	for (int i = 0; i != fnum; ++i)
	{
		fixMat[i](0) = fix[i][0];
		fixMat[i](1) = fix[i][1];
		fixMat[i](2) = fix[i][2];
	}
	char strPara[MAX_PATH] = { 0 };
	int len = GetPrivateProfileStringA("Param",
		"Deg",
		"No Deg",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	int deg = atoi(strPara);

	len = GetPrivateProfileStringA("Param",
		"maxIterCount",
		"No maxIterCount",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	int maxIterCount = atoi(strPara);

	len = GetPrivateProfileStringA("Param",
		"Omega",
		"No Omega",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	double omega = atof(strPara);

	len = GetPrivateProfileStringA("Param",
		"Eps",
		"No Eps",
		strPara,
		MAX_PATH,
		".//ACPD.ini");

	double eps = atof(strPara);

	char error_time[MAX_PATH] = { 0 };
	len = GetPrivateProfileStringA("Param",
		"OutputErrorTime",
		"No OutputErrorTime",
		error_time,
		MAX_PATH,
		".//ACPD.ini");

	double error = 0;
	error = Analytic_CPD_3D(
		moMat,
		mnum,
		fixMat,
		fnum,
		maxIterCount,
		eps,
		deg,
		omega,
		error_time);

	for (int i = 0; i != mnum; ++i)
	{
		mo[i][0] = moMat[i](0);
		mo[i][1] = moMat[i](1);
		mo[i][2] = moMat[i](2);
	}

	delete[] moMat;
	delete[] fixMat;

	return error;
}

DLL_ANALYTICCPD double feature_transf2d(double(*mov)[2], int mnum
	, double(*fix)[2], int fnum, int iterCount)
{
	Eigen::setNbThreads(1);

	std::cout << "Eigen threads = "
		<< Eigen::nbThreads()
		<< std::endl;
	double error = feature_nonrigid_regist2d(mov, mnum, fix, fnum, iterCount);

	return error;
}

DLL_ANALYTICCPD double feature_transf3d(double(*mo)[3], int mnum
	, double(*fix)[3], int fnum, int iterCount)
{
	Eigen::setNbThreads(1);

	std::cout << "Eigen threads = "
		<< Eigen::nbThreads()
		<< std::endl;
	double error = feature_nonrigid_regist3d(mo, mnum, fix, fnum, iterCount);

	return error;
}

DLL_ANALYTICCPD void randomly_2d_map(double(*taget)[2], double(*moved)[2]
	, int num, int deg)
{
	AnalyticTransf2d(taget, moved, num, deg);
}

DLL_ANALYTICCPD void randomly_3d_map(double(*target)[3], double(*moved)[3]
	, int num, int deg)
{
	AnalyticT3dArrayWrapper(target, moved, num, deg, static_cast<unsigned int>(time(nullptr)));
}