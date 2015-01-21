#include "timer.hpp"
#include <stdint.h>


timestamp_t getTimeinMilliSeconds(){
    struct timeval timer;
    gettimeofday(&timer, NULL);
    //Number of microseconds

    timestamp_t returnTime = timer.tv_usec;

    //Number of milliseconds;
    returnTime /= 1000.0;

    //Total millisecondsl
    returnTime += (timer.tv_sec*1000.0);
    return returnTime;

}
