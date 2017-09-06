#include <stdio.h>
#include <iostream>
#include <cstdint>
#include "triu_index.h"

int main (int argc, char *argv[])
{

	std::cerr << false << std::endl;
	Triangular_index T(4);
	size_t N=4;
	
	for (size_t x=0; x<N; ++x)
	{
		for (size_t y=x+1; y<N; ++y)
		{
			for (size_t z=0; z<N; ++z)
			{
				for (size_t a=0; a<3; ++a)
				{
					for (size_t b=0; b<3; ++b)
					{
						std::cerr << x << ", " << y << ", " << z << ", " << a << ", " << b << "->" << T.get_k(x,y)*9*N+9*z+a*3+b << std::endl;
					}
				}
				
			}
		}
	}
	for (Triangular_index T(4); T.get_k()<T.size(); T++)
	{
		size_t x, y;
		T.get_xy(x,y);
		fprintf(stderr, "%ld, %ld, %ld\n", x, y, T.get_k() );
	}
}
