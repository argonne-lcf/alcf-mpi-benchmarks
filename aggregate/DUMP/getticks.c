
#include <hwi/include/bqc/A2_inlines.h>

#define SEC_PER_CYCLE                                      (1.0 / 1600000000.0)
double timer()
{
  return ((double)GetTimeBase() * SEC_PER_CYCLE);
}

unsigned long long getticks(void)
{
     unsigned int rx, ry, rz;
     unsigned long long r64;

     do
     {
         asm volatile ( "mftbu %0" : "=r"(rx) );
         asm volatile ( "mftb %0" : "=r"(ry) );
         asm volatile ( "mftbu %0" : "=r"(rz) );
     }
     while ( rx != rz );

     r64 = rx;
     r64 = ( r64 << 32 ) | ry;
     
     return r64;
}

