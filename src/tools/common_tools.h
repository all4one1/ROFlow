#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <windows.h>

using std::string;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::cout;
using std::endl;


struct ReadingFile
{
private:
	ifstream read;
	string str, substr, buffer;
	string file_name;
	stringstream ss;
	istringstream iss;
	ostringstream oss;
	int stat = 0, pos = 0;

public:
	ReadingFile(string name)
	{
		file_name = name;
		open_file(file_name);
		stat = 0;
	}
	ReadingFile()
	{
		stat = 0;
	}

	void open_file(string file_name) {
		read.open(file_name.c_str());
		if (read.good()) {
			cout << endl << "the parameter file \"" << file_name << "\" has been read " << endl << endl;
			#ifdef LoggingON
			Log << endl << "the parameter file \"" << file_name << "\" has been read " << endl << endl;
			#endif // LoggingON

			oss << read.rdbuf();
			buffer = oss.str();
			iss.str(buffer);
		}
		else {
			cout << "the parameter file has been not found, default parameters will be initialized " << endl;
			#ifdef LoggingON
			Log << "the parameter file has been not found, default parameters will be initialized " << endl;
			#endif // LoggingON
			buffer = "";
			iss.str(buffer);
		}
	}

	template <typename T>
	int reading(std::map<string, T>& m, T def_var = 0, T min = 0, T max = 0)
	{
		for (auto& it : m)
		{
			reading<T>(it.second, it.first, def_var, min, max);
		}
		return 0;
	}
	template <typename T>
	int reading(T& var, string parameter_name, T def_var = 0, T min = 0, T max = 0) {
		int ret = 0;
		stat = 0;
		transform(parameter_name.begin(), parameter_name.end(), parameter_name.begin(), ::tolower);
		iss.clear();
		iss.seekg(0);

		while (getline(iss, str))
		{
			//substr.clear();
			ss.str("");	ss.clear();	ss << str;	ss >> substr;
			transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
			if (substr == parameter_name) {
				ret = 1;
				pos = (int)ss.tellg();
				while (ss >> substr) {
					if (substr == "=")
					{
						ss >> var;
						stat = 1;
						break;
					}
				}

				if (stat == 0) {
					ss.clear();
					ss.seekg(pos);
					ss >> var;
				}
				break;
			}
		}
		if (iss.fail())
		{
			var = def_var;
		}
		if (min != max && (min + max) != 0) {
			if (var > max || var < min)
			{
				cout << "Warning: \"" + parameter_name + "\" should not be within this range" << endl;
				var = def_var;
			}
		}


		#ifdef LoggingON
		Log << parameter_name << "= " << var << endl;
		#endif // LoggingON

		return ret; //return 1 if read
	}

	void reading_string(string& var, string parameter_name, string def_var) {
		stat = 0;
		transform(parameter_name.begin(), parameter_name.end(), parameter_name.begin(), ::tolower);
		iss.clear();
		iss.seekg(0);

		while (getline(iss, str))
		{
			//substr.clear();
			ss.str("");	ss.clear();	ss << str;	ss >> substr;
			transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
			if (substr == parameter_name) {
				pos = (int)ss.tellg();
				while (ss >> substr) {
					if (substr == "=")
					{
						ss >> var;
						stat = 1;
						break;
					}
				}

				if (stat == 0) {
					ss.clear();
					ss.seekg(pos);
					ss >> var;
				}
				break;
			}
		}
		if (iss.fail())
		{
			var = def_var;
		}
		#ifdef LoggingON
		Log << parameter_name << "= " << var << endl;
		#endif // LoggingON
	}
	bool reading_string(string parameter_name) {
		stat = 0;
		transform(parameter_name.begin(), parameter_name.end(), parameter_name.begin(), ::tolower);
		iss.clear();
		iss.seekg(0);

		while (getline(iss, str))
		{
			//substr.clear();
			ss.str("");	ss.clear();	ss << str;	ss >> substr;
			transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
			if (substr == parameter_name)
			{
				return true;
			}
		}
		#ifdef LoggingON
		Log << parameter_name << "= " << var << endl;
		#endif // LoggingON
		return false;
	}
};

struct FuncTimer
{
private:
	std::map <std::string, double> timer, t1, t2;
	std::string active;
	double dt, total;
public:
	FuncTimer()
	{
		dt = total = 0;
	}
	void start(std::string s)
	{
		active = s;
		t1[s] = clock();
	}
	void end(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return;
		}
		else
		{
			t2[s] = clock();
			dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
			//if (s == active)
			{
				timer[s] += dt;
			}
		}
	}
	double get(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return 0.0;
		}
		else
		{
			return timer[s];
		}
	}

	double update_and_get(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return 0.0;
		}
		else
		{
			t2[s] = clock();
			dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
			return timer[s] + dt;
		}
	}
	double get_last_diff(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return 0.0;
		}
		else
		{
			t2[s] = clock();
			dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
			return dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
		}
	}

	std::string get_info()
	{
		int n = int(timer.size());
		std::ostringstream oss;
		oss << "Calculation time in seconds. Number of cases: " << n << ".\n";

		for (auto& it : timer)
		{
			oss << it.first << ": " << it.second << std::endl;
		}
		return oss.str();
	}
	void show_info()
	{
		std::cout << get_info() << std::endl;
	}
	void write_info(std::string path = "timer_report.dat")
	{
		std::ofstream file(path);
		file << get_info() << std::endl;
	}

};

struct Checker
{
	enum class ExitType { Relative, Absolute } type;
	std::string name;
	std::vector<double> v, t;
	double* binded_value = nullptr, * binded_time = nullptr;
	double t_estimated = 0.0;
	double last_time = 0;
	size_t N = 0, i2 = 0, i1 = 0, i0 = 0, iter = 0;
	double eps_default = 1e-6;
	double relative = 1, absolute = 1;
	double dif = 1, dif_rel = 1, dif2 = 0, dif_rel2 = 0;
	double deriv = 0, deriv2 = 0;
	size_t i_check = 0, check_limit = 10;
	bool ready_to_exit = false;

	Checker(double* ptr, double* t_ptr, ExitType type_ = ExitType::Relative, std::string name_ = "Checker", double eps = -1)
	{
		iter = 0;
		N = 3;
		binded_value = ptr;
		binded_time = t_ptr;
		v.resize(N);
		std::fill(v.begin(), v.end(), 0);

		t.resize(N);
		std::fill(t.begin(), t.end(), 0);

		type = type_;
		name = name_;

		if (eps > 0) eps_default = eps;
	}

	void update()
	{
		if (last_time == binded_time[0]) return;
		i2 = iter % N;
		i1 = (iter - 1) % N;
		i0 = (iter - 2) % N;

		v[i2] = binded_value[0];
		t[i2] = binded_time[0];
		iter++;
		last_time = binded_time[0];
		if (iter < N) return;

		dif = v[i2] - v[i1];
		dif_rel = (v[i2] - v[i1]) / v[i1];
		//dif2 = (v[i2] - v[i1]) - (v[i1] - v[i0]);
		//dif_rel2 = ((v[i2] - v[i1]) - (v[i1] - v[i0])) / (v[i2] - v[i1]);
		double dt = t[i2] - t[i1];
		deriv = (v[i2] - v[i1]) / dt;
		deriv2 = (v[i2] - 2.0 * v[i1] + v[i0]) / (dt * dt);


		absolute = abs(dif);
		relative = abs(dif / v[i1]);

		check();
		//cout << iter << ": " << i2 << " " << i1 << " " << i0 << endl;
		//cout << iter << ": " << v[i2] << " " << v[i1] << " " << v[i0] << endl;
	}

	void update(double& val, double& t)
	{
		binded_value = &val;
		binded_time = &t;
		this->update();
	}

	bool check(double eps = -1)
	{
		//update();
		if (eps < 0) eps = eps_default;
		if (relative < eps && type == ExitType::Relative
			|| absolute < eps && type == ExitType::Absolute)
		{
			i_check++;
			cout << name + ",i= " << i_check << endl;
			if (i_check >= check_limit)
			{
				ready_to_exit = true;
				return true;
			}
		}
		else
		{
			i_check = 0;
		}
		ready_to_exit = false;
		return false;
	}

	double estimate_finish(double eps = -1)
	{
		update();
		if (eps < 0) eps = eps_default;

		double e2 = (v[i2] - v[i1]) / v[i1];
		double e1 = (v[i1] - v[i0]) / v[i0];
		double dt = t[i2] - t[i1];

		return t[i2] + dt * (eps - e2) / (e2 - e1);
	}

	bool operator()()
	{
		return ready_to_exit;
	}

	//bool check_long_term();
};


inline size_t allocate_host_arrays(std::vector<double**> v, size_t N)
{
	size_t Nbytes = sizeof(double) * N;
	size_t total_bytes = 0;
	for (auto it : v)
	{
		*it = new double[N]();
		//memset(*it, 0, Nbytes);
		total_bytes += Nbytes;
	}
	return total_bytes;
}
inline void addFolder(std::map<string, string>& m, std::vector<string> names)
{
	for (auto& it : names)
	{
		#ifdef __linux__ 
		string str = "mkdir -p " + it + "/";
		m[it] = it + "/";
		system(str.c_str());
		#endif

		#ifdef _WIN32
		CreateDirectoryA(it.c_str(), NULL);
		m[it] = it + "//";
		#endif
	}
}
inline void deleteFilesInDirectory(const std::wstring& directoryPath) {
	WIN32_FIND_DATAW findFileData;
	HANDLE hFind = INVALID_HANDLE_VALUE;

	// Create a search path for all files in the directory
	std::wstring searchPath = directoryPath + L"\\*";

	// Find the first file in the directory
	hFind = FindFirstFileW(searchPath.c_str(), &findFileData);
	if (hFind == INVALID_HANDLE_VALUE) {
		std::wcerr << L"Error opening directory: " << directoryPath << std::endl;
		return;
	}

	do {
		// Skip the special directories "." and ".."
		if (wcscmp(findFileData.cFileName, L".") == 0 || wcscmp(findFileData.cFileName, L"..") == 0) {
			continue;
		}

		// Construct the full file path
		std::wstring filePath = directoryPath + L"\\" + findFileData.cFileName;

		// Check if the entry is a regular file
		if (!(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
			// Delete the file
			if (DeleteFileW(filePath.c_str()))
			{
				//std::wcout << L"Deleted file: " << filePath << std::endl;
			}
			else {
				std::wcerr << L"Error deleting file: " << filePath << std::endl;
			}
		}
	} while (FindNextFileW(hFind, &findFileData) != 0);

	// Close the search handle
	FindClose(hFind);
}
inline std::streamoff fileSize(const std::string& filePath) {
	std::ifstream file(filePath, std::ios::binary | std::ios::ate);
	if (!file.is_open()) { std::cerr << "Error!" << std::endl;        return false; }
	std::streampos bytes = file.tellg();
	return bytes.operator std::streamoff(); 	//file size in bytes
}


inline double sum_signed(double* f, unsigned int N)
{
	double s = 0;

	for (unsigned int i = 0; i < N; i++)
		s += f[i];
	return s;
}
inline double sum_abs(double* f, unsigned int N)
{
	double s = 0;

	for (unsigned int i = 0; i < N; i++)
		s += abs(f[i]);
	return s;
}
inline double max_signed(double* f, unsigned int N)
{
	double m = f[0];

	for (unsigned int i = 0; i < N; i++)
	{
		if (f[i] > m) m = f[i];
	}
	return m;
}
inline double min_signed(double* f, unsigned int N)
{
	double m = f[0];

	for (unsigned int i = 0; i < N; i++)
	{
		if (f[i] < m) m = f[i];
	}
	return m;
}
inline double max_abs(double* f, unsigned int N)
{
	double m = f[0];

	for (unsigned int i = 0; i < N; i++)
	{
		if (abs(f[i]) > m) m = abs(f[i]);
	}
	return m;
}

inline double minval_signed(const std::vector<double>& val)
{
	if (val.empty()) {
		throw std::invalid_argument("Vector is empty");
	}
	return *std::min_element(val.begin(), val.end());
}
inline double minval_abs(const std::vector<double>& val) {
	if (val.empty()) {
		throw std::invalid_argument("Vector is empty");
	}
	auto it = std::min_element(val.begin(), val.end(),
		[](double a, double b) {
		return std::abs(a) < std::abs(b);
	});
	return *it;
}

inline std::string _str(double n, int pres = 5)
{
	std::stringstream ss;
	ss << std::setprecision(pres) << n;
	return ss.str();
}
inline std::string _str(int n, int pres = 5)
{
	std::stringstream ss;
	ss << std::setprecision(pres) << n;
	return ss.str();
}
inline std::string _path(std::string fname)
{
	#ifdef __linux__ 
	return fname + "/";
	#endif

	#ifdef _WIN32
	return fname + "//";
	#endif
}
