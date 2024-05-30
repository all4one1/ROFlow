#define __FILENAME__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#define print(message) {std::cout << __FILENAME__  << " (Line #" << __LINE__ << "): " << message << std::endl;}
#define pause {std::cout << "Pause at line: " << __LINE__ << " in file: " << __FILENAME__  << std::endl; system("pause"); }
#define MYERROR(message) print(message); throw std::invalid_argument(message); 


void allocM(double** f, int N);
void alloc(double** f, int N);