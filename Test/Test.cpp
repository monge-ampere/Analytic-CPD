// test.cpp : Console frontend for Analytic-CPD
//

#include "stdafx.h"
#include <iostream>
#include <windows.h>
#include <tchar.h>
#include <cstdlib>

using namespace std;

static bool ReadIniString(
	const char* section,
	const char* key,
	const char* defaultValue,
	char* output,
	DWORD outputSize,
	const char* iniPath)
{
	DWORD len = GetPrivateProfileStringA(
		section,
		key,
		defaultValue,
		output,
		outputSize,
		iniPath);

	return len > 0;
}

static int ReadIniInt(
	const char* section,
	const char* key,
	int defaultValue,
	const char* iniPath)
{
	char buf[MAX_PATH] = { 0 };
	char defaultBuf[64] = { 0 };
	sprintf_s(defaultBuf, "%d", defaultValue);

	GetPrivateProfileStringA(
		section,
		key,
		defaultBuf,
		buf,
		MAX_PATH,
		iniPath);

	return atoi(buf);
}

int main()
{
	typedef bool(*ma_initialize)(TCHAR* dllPath);

	typedef double(*ma_ps_regist)(
		const char* fixed_path,
		const char* moving_path,
		const char* moved_output);

	typedef void(*ma_analytic_map)(
		const char* fixed_path,
		const char* deformed_output,
		int deg);

	typedef void(*ma_release)();

	const char* experimentIni = ".\\experiment.ini";

	printf("Analytic-CPD console frontend begins ----------------------------\n\n");

	// ---------------------------------------------------------------------
	// Load SmoothAdjustment.dll
	// ---------------------------------------------------------------------
	HINSTANCE hDll = LoadLibrary(L"SmoothAdjustment.dll");
	if (!hDll)
	{
		cerr << "Error: failed to load SmoothAdjustment.dll." << endl;
		cerr << "Please make sure SmoothAdjustment.dll is in the runtime directory." << endl;
		system("pause");
		return -1;
	}

	ma_initialize initialize =
		(ma_initialize)GetProcAddress(hDll, "acpd_initialize");

	ma_release release =
		(ma_release)GetProcAddress(hDll, "acpd_release");

	ma_ps_regist regist2d =
		(ma_ps_regist)GetProcAddress(hDll, "acpd_register_2d_from_file");

	ma_ps_regist regist3d =
		(ma_ps_regist)GetProcAddress(hDll, "acpd_register_3d_from_file");

	ma_analytic_map analyticMap2d =
		(ma_analytic_map)GetProcAddress(hDll, "analytic_2d_map");

	ma_analytic_map analyticMap3d =
		(ma_analytic_map)GetProcAddress(hDll, "analytic_3d_map");

	if (!initialize || !release || !regist2d || !regist3d)
	{
		cerr << "Error: failed to load required functions from SmoothAdjustment.dll." << endl;
		FreeLibrary(hDll);
		system("pause");
		return -1;
	}

	// ---------------------------------------------------------------------
	// Read experiment.ini
	// ---------------------------------------------------------------------
	int registType = ReadIniInt("Param", "RegistType", 1, experimentIni);
	int addPerturb = ReadIniInt("Param", "AddPerturb", 0, experimentIni);
	int perturbDeg2D = ReadIniInt("Param", "PerturbDeg2D", 8, experimentIni);
	int perturbDeg3D = ReadIniInt("Param", "PerturbDeg3D", 2, experimentIni);

	char movingPsPath[MAX_PATH] = { 0 };
	char fixedPsPath[MAX_PATH] = { 0 };
	char movedOutputPath[MAX_PATH] = { 0 };
	char deformedOutputPath[MAX_PATH] = { 0 };

	ReadIniString("Param", "MovingPsPath", "", movingPsPath, MAX_PATH, experimentIni);
	ReadIniString("Param", "FixedPsPath", "", fixedPsPath, MAX_PATH, experimentIni);
	ReadIniString("Param", "MovedOutputPath", "results\\moved.csv", movedOutputPath, MAX_PATH, experimentIni);
	ReadIniString("Param", "DeformedOutputPath", "results\\deformed.csv", deformedOutputPath, MAX_PATH, experimentIni);

	if (strlen(movingPsPath) == 0 || strlen(fixedPsPath) == 0)
	{
		cerr << "Error: MovingPsPath or FixedPsPath is empty in experiment.ini." << endl;
		FreeLibrary(hDll);
		system("pause");
		return -1;
	}

	// ---------------------------------------------------------------------
	// Initialize Analytic_CPD.dll through the wrapper
	// ---------------------------------------------------------------------
	if (!initialize(_T("Analytic_CPD.dll")))
	{
		cerr << "Error: failed to initialize Analytic_CPD.dll." << endl;
		FreeLibrary(hDll);
		system("pause");
		return -1;
	}

	printf("Non-rigid registration with Analytic-CPD begins -----------------\n");
	printf("Fixed point set : %s\n", fixedPsPath);
	printf("Moving point set: %s\n", movingPsPath);
	printf("Output point set: %s\n", movedOutputPath);

	double result = -1.0;

	if (registType == 1)
	{
		// 3D registration
		if (addPerturb)
		{
			if (!analyticMap3d)
			{
				cerr << "Error: analytic_3d_map is not available." << endl;
				release();
				FreeLibrary(hDll);
				system("pause");
				return -1;
			}

			printf("Generating 3D analytic perturbation: degree = %d\n", perturbDeg3D);
			analyticMap3d(movingPsPath, deformedOutputPath, perturbDeg3D);
		}

		result = regist3d(fixedPsPath, movingPsPath, movedOutputPath);
	}
	else
	{
		// 2D registration
		if (addPerturb)
		{
			if (!analyticMap2d)
			{
				cerr << "Error: analytic_2d_map is not available." << endl;
				release();
				FreeLibrary(hDll);
				system("pause");
				return -1;
			}

			printf("Generating 2D analytic perturbation: degree = %d\n", perturbDeg2D);
			analyticMap2d(movingPsPath, deformedOutputPath, perturbDeg2D);
		}

		result = regist2d(fixedPsPath, movingPsPath, movedOutputPath);
	}

	printf("\nRegistration finished.\n");
	printf("Final error: %.10f\n", result);

	release();
	FreeLibrary(hDll);

	printf("Analytic-CPD console frontend finishes --------------------------\n");
	system("pause");
	return 0;
}