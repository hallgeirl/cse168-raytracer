#ifndef WIN32
#include <sys/time.h>
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
    printf("getTime() not implemented for this platform.\n");
    //TODO: Implement a Windows version of this.
    return 0;
    #endif
}

