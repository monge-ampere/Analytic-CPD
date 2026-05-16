// Analytic_CPD.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "Analytic_CPD.h"


// 这是导出变量的一个示例
ANALYTIC_CPD_API int nAnalytic_CPD=0;

// 这是导出函数的一个示例。
ANALYTIC_CPD_API int fnAnalytic_CPD(void)
{
	return 42;
}

// 这是已导出类的构造函数。
// 有关类定义的信息，请参阅 Analytic_CPD.h
CAnalytic_CPD::CAnalytic_CPD()
{
	return;
}
