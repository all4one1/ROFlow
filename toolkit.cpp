#pragma once
#include <malloc.h>
#include <iostream>

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
void alloc(double** f, int N)
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
};