// swapping stringstream objects
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include "state.h"

int main () {

	State a(1);
	State b(3);
	uint32_t m[1]={0x5};
	uint32_t P1[5]={0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555}, P2[5]={0x66666666, 0x77777777, 0x88888888, 0x99999999, 0xAAAAAAAA};
	for (size_t x=0; x<5; x++) a.compress(P1+x, P1+x);
	a.cache();
	a.rewind();
	for (size_t x=0; x<5; x++) 
	{	
		a.uncompress(P1, P2);
		std::cerr << P1[0] << std::endl;
	}
	std::cerr << std::endl;
	a.rewind();
	for (size_t x=5; x>0; x--) 
	{	
		a.uncompress(P1, P2, x);
		std::cerr << P1[0] << std::endl;
	}
	return 0;
}
