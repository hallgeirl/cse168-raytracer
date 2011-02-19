#ifndef WIN32
#include <sys/time.h>
#else
#include <windows.h>
#endif
#include <stdio.h>

//Returns the time in seconds
double getTime()
{
    #ifndef WIN32
    struct timeval tv;
    int ret = gettimeofday(&tv, 0);
    double t = (double)tv.tv_sec + 1e-6 * (double)tv.tv_usec;
    if (ret == -1)
    {
        printf("Something bad happened when calling gettimeofday().\n");
        return -1;
    }
    return t;

    #else
	return GetTickCount()/1000.f;
	#endif
}

