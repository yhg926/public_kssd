#ifndef MYTIME
#define MYTIME 
long long timeInMilliseconds(void);
long long timeInMicroseconds(void);
struct timespec diff(struct timespec start, struct timespec end);
#endif
