// 下列 ifdef 块是创建使从 DLL 导出更简单的
// 宏的标准方法。此 DLL 中的所有文件都是用命令行上定义的 ANALYTIC_CPD_EXPORTS
// 符号编译的。在使用此 DLL 的
// 任何其他项目上不应定义此符号。这样，源文件中包含此文件的任何其他项目都会将
// ANALYTIC_CPD_API 函数视为是从 DLL 导入的，而此 DLL 则将用此宏定义的
// 符号视为是被导出的。
#ifdef ANALYTIC_CPD_EXPORTS
#define ANALYTIC_CPD_API __declspec(dllexport)
#else
#define ANALYTIC_CPD_API __declspec(dllimport)
#endif

// 此类是从 Analytic_CPD.dll 导出的
class ANALYTIC_CPD_API CAnalytic_CPD {
public:
	CAnalytic_CPD(void);
	// TODO:  在此添加您的方法。
};

extern ANALYTIC_CPD_API int nAnalytic_CPD;

ANALYTIC_CPD_API int fnAnalytic_CPD(void);
