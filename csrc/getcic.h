#ifndef POPHIST_GETCIC_H
#define POPHIST_GETCIC_H

#define START_COMMENT '#'

int getcic(FILE * fp);
unsigned getucic(FILE * fp);
char *getwordic(char *buff, int bufsiz, FILE * ifp);
int getdblic(double *x, char *buff, int bufsiz, FILE * ifp);
int getintic(int *i, char *buff, int bufsiz, FILE * ifp);
int getunsignedic(unsigned *u, char *buff, int bufsiz, FILE * ifp);

#endif
