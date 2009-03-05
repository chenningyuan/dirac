#include <time.h>

clock_t mclock_(void)
{
clock_t tmp;

tmp = clock();
return tmp;
}
