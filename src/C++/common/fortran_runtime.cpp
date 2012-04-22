#include <cstdio>
// These appear to be needed in order to disable linker errors about
// missing Fortran runtime functions
extern "C" void for_write_seq_fmt(...)
{
    printf("%s\n", __FUNCTION__); fflush(stdout);
}

extern "C" void for_write_seq_fmt_xmit(...)
{
    printf("%s\n", __FUNCTION__); fflush(stdout);
}