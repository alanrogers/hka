
/*
 * The functions handle errors and then terminate execution.  They are copied
 * from pp 109-114 of:
 *  Brian Kernighan and Rob Pike. 1997. The Practice of Programming. 
 *  Addison-Wesley.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include "eprintf.h"

/* eprintf: print error message and exit */
void
eprintf(const char *fmt, ...)
{
    va_list args;

    fflush(stdout);

    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    va_end(args);

    if(fmt[0] != '\0' && fmt[strlen(fmt) - 1] == ':')
        fprintf(stdout, " %s", strerror(errno));
    fprintf(stdout, "\n");
    exit(2);                    /* conventional value for failed execution */
}

#ifdef XEPRINTF
#include <unistd.h>
#include <math.h>

/* test the eprintf function */
void
main(void)
{
    int i = 2;
    double x = 3.4;

    if(fork() == 0) {
        sleep(5);
        eprintf("eprintf w/ no args:");
    }

    if(fork() == 0) {
        eprintf("eprintf w/ no args:");
    }

    if(fork() == 0) {
        eprintf("eprintf w/ args. i=%d x=%lf:", i, x);
    }

    if(fork() == 0) {
        x /= 0.0;
        eprintf("eprintf after division by 0. x=%lf:", x);
    }

    if(fork() == 0) {
        x = log(0.0);
        eprintf("eprintf after log(0.0) x=%lf:", x);
    }

    sleep(1);

    printf("Parent exiting with status 0\n");
    exit(0);
}
#endif
