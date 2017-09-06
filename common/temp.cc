#include <unistd.h>

int main (int argc, char *argv [])
{
	int *c=new int [1000000000]{11};
	usleep(10000000);
	return 0;
}
