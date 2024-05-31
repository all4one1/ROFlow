#pragma once
#include <malloc.h>
#include <iostream>
#include <vector>
#pragma warning( disable: 6386 )
void allocM(double** f, int N)
{
	*f = (double*)malloc(sizeof(double) * N);
	if (*f == nullptr)
	{
		std::cout << "memory allocation problem!" << std::endl;
	}
	else
	{
		for (int l = 0; l < N; l++)
		{
			(*f)[l] = 0.0;
		}
	}
};
size_t alloc(double** f, int N)
{
	try
	{
		*f = new double[N];
		for (int l = 0; l < N; l++)
		{
			(*f)[l] = 0.0;
		}
	}
	catch (const std::bad_alloc& e)
	{
		std::cout << "Allocation failed: " << e.what() << '\n';
	}

	return sizeof(double) * N;
};

struct MemoryAllocator
{
	size_t bytes_ = 0;

	void alloc(double** f, int N)
	{
		try
		{
			*f = new double[N];
			for (int l = 0; l < N; l++)
			{
				(*f)[l] = 0.0;
			}
			bytes_ += sizeof(double) * N;
		}
		catch (const std::bad_alloc& e)
		{
			std::cout << "Allocation failed: " << e.what() << '\n';
		}
	};

	void allocvec(int N, std::vector<double**> v)
	{
		for (auto& it : v)
		{

			*it = new double[N];
			for (int l = 0; l < N; l++)
			{
				(*it)[l] = 0.0;
			}

		}
	}

};