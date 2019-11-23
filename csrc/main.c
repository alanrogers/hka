/**
   @file hka.c
   @brief Perform HKA test of selective neutrality
   @internal
   Copyright (C) 2003 Alan R. Rogers

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. rogers@anthro.utah.edu
**/
#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  "hka.h"
#include  "getcic.h"
#include  "eprintf.h"
#include  "hkaversion.h"

void    usage(void);
char  * strndup(const char *s, size_t n);

void    usage(void)
{
    printf("usage: hka filename\n");
    exit(1);
}

/*
 * Like strdup, but duplicated string contains at most n characters
 * including the terminating '\0'.
 */
char  * strndup(const char *s, size_t n) {
    size_t len;
    char * u;

    len = strlen(s);
    if(len > n-1)
	len = n-1;

    u = malloc(len+1);
    memcpy(u, s, len);
    u[len] = '\0';
    return u;
}

int main(int argc, char **argv)
{
    FILE *fp = NULL;
    unsigned i, nLoci;
    int buffsize = 100;
    char * buff;

    char ** locusName;
    unsigned * nsA;   /* number of sites A */
    unsigned * nsB;   /* number of sites B */
    unsigned * nsAB;  /* number of sites AB */
    unsigned * nA;    /* haploid sample size A */
    unsigned * nB;    /* haploid sample size B */
    double   * SA;    /* polymorphic sites A */
    double   * SB;    /* polymorphic sites B */
    double   * D;     /* mean diff AB */
    HKA * h;          /* object that does the analysis */

    buff = (char *) malloc((unsigned) buffsize);
    if(buff == NULL)
	eprintf("main: no memory\n");

    /* Process command line arguments */
    if(argc > 1) {
	fp = fopen(argv[1], "r");
    }

    if(fp == NULL)
	usage();

    // Read number of loci
    getunsignedic(&nLoci, buff, buffsize, fp);

    // Allocate arrays
    locusName = (char **) malloc( nLoci * sizeof(char *) );
    nsA       = (unsigned *) malloc(nLoci * sizeof( unsigned ) );
    nsB       = (unsigned *) malloc(nLoci * sizeof( unsigned ) );
    nsAB      = (unsigned *) malloc(nLoci * sizeof( unsigned ) );
    nA        = (unsigned *) malloc(nLoci * sizeof( unsigned ) );
    nB        = (unsigned *) malloc(nLoci * sizeof( unsigned ) );
    SA        = (double *) malloc(nLoci * sizeof( double ) );
    SB        = (double *) malloc(nLoci * sizeof( double ) );
    D         = (double *) malloc(nLoci * sizeof( double ) );

    if(locusName==NULL || nsA==NULL || nsB==NULL || nA==NULL || nB==NULL
       || SA==NULL || SB==NULL || D==NULL)
	eprintf("main: no memory\n");

    /* Finish reading input */
    for(i=0; i < nLoci; ++i) {
	getwordic(buff, buffsize, fp);
	locusName[i] = strndup(buff, (size_t) buffsize);
	getunsignedic(nsA+i, buff, buffsize, fp);
	getunsignedic(nsB+i, buff, buffsize, fp);
	getunsignedic(nsAB+i, buff, buffsize, fp);
	getunsignedic(nA+i, buff, buffsize, fp);
	getunsignedic(nB+i, buff, buffsize, fp);
	getdblic(SA+i, buff, buffsize, fp);
	getdblic(SB+i, buff, buffsize, fp);
	getdblic(D+i, buff, buffsize, fp);
    }

    printf("\n########## HKA Version %s  ##########\n",
	   HKA_VERSION);

    /* Echo command line */
    printf("# Command line:");
    for(i=0; i < argc; ++i)
	printf(" %s", argv[i]);
    putchar('\n');

    /* Echo data */
    printf("%d # loci\n", nLoci);
    printf("%10s %6s %6s %6s %6s %6s %6s %6s %6s\n",
	   "#    Locus", "nsA", "nsB", "nsAB", "nA", "nB",
	   "SA", "SB", "D");

    for(i=0; i < nLoci; ++i) {
	printf("%10s %6u %6u %6u %6u %6u %6g %6g %6g\n",
	       locusName[i], nsA[i], nsB[i], nsAB[i],
	       nA[i], nB[i], SA[i], SB[i], D[i]);
    }
    puts("#");
    puts("# nsA  Number of sites used w/i species A");
    puts("# nsB  w/i species B");
    puts("# nsAB in comparison btw A and B");
    puts("# nA   Haploid sample size for A");
    puts("# nB   for B");
    puts("# SA   Polymorphic sites w/i A");
    puts("# SB   w/i B");
    puts("# D    Mean number of site differences btw A and B");

    /* Create object of type HKA */
    h = HKA_construct(nLoci, nA, nB, nsA, nsB, nsAB, SA, SB, D);

    /* Print results */
    puts("\n## RESULTS ##");
    printf("#        T = %f\n", HKA_T(h));
    printf("#        f = %f\n", HKA_f(h));
    for(i=0; i < nLoci; ++i)
	printf("# theta[%d] = %f\n", i, HKA_theta(i, h));
    printf("#       df = %u\n", HKA_degreesOfFreedom(h));
    printf("#      X^2 = %f\n", HKA_X2(h));

    return 0;
}
