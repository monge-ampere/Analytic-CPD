#pragma once
#include <fstream>
#include <iomanip>
#include <string>
#include <stdexcept>

struct RegistrationIterLog
{
	int iter = 0;
	int degree = 0;

	double rmse = 0.0;       
	double meanL2 = 0.0;      
	double softRMSE = 0.0;    
	double sigma2 = 0.0;
	double fitErr = 0.0;

	double runtimeSec = 0.0;

	int activeCount = 0;
	double NP = 0.0;
};

class RegistrationCSVLogger
{
public:
	RegistrationCSVLogger() = default;

	explicit RegistrationCSVLogger(
		const std::string& filename,
		bool detailed = false,
		bool header = false)
	{
		open(filename, detailed, header);
	}

	void open(
		const std::string& filename,
		bool detailed = false,
		bool header = false)
	{
		detailed_ = detailed;
		fp_.open(filename.c_str(), std::ios::out);

		if (!fp_.is_open()) {
			throw std::runtime_error("RegistrationCSVLogger: failed to open file: " + filename);
		}

		fp_ << std::setprecision(15);

		if (header) {
			if (detailed_) {
				fp_ << "iter,degree,rmse,meanL2,softRMSE,sigma2,fitErr,runtimeSec,activeCount,NP\n";
			}
			else {
				fp_ << "rmse,runtimeSec\n";
			}
		}
	}

	bool is_open() const
	{
		return fp_.is_open();
	}

	void write(const RegistrationIterLog& r)
	{
		if (!fp_.is_open()) {
			return;
		}

		if (detailed_) {
			fp_
				<< r.iter << ","
				<< r.degree << ","
				<< r.rmse << ","
				<< r.meanL2 << ","
				<< r.softRMSE << ","
				<< r.sigma2 << ","
				<< r.fitErr << ","
				<< r.runtimeSec << ","
				<< r.activeCount << ","
				<< r.NP
				<< "\n";
		}
		else {
			fp_ << r.rmse << "," << r.runtimeSec << "\n";
		}
	}

	void close()
	{
		if (fp_.is_open()) {
			fp_.close();
		}
	}

	~RegistrationCSVLogger()
	{
		close();
	}

private:
	std::ofstream fp_;
	bool detailed_;
};