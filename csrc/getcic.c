#include <stdio.h>
#include <ctype.h>
#include "getcic.h"

/*
 * Get an unsigned char from file stream fp, ignoring comments.
 */
unsigned
getucic(FILE * fp)
{
    return (unsigned) getcic(fp);
}

/*
 * Modify getc() to ignore comments, which begin with START_COMMENT and end
 * with newline or EOF.  When a comment appears getcic() returns the
 * newline (or EOF) character only.
 */
int
getcic(FILE * fp)
{
    int c;

    c = getc(fp);
    if(c == START_COMMENT) {
        do {
            c = getc(fp);       /* read to end of line */
        } while(c != '\n' && c != EOF);
    }
    return (c);
}

/*
 * get a word, delimited by whitespace, from file ifp
 * Ignore comments delimited by START_COMMENT...\n
 */
char *
getwordic(char *buff, int bufsiz, FILE * ifp)
{
    char *bp;
    int c;

    bp = buff;
    do {
        c = getcic(ifp);
    } while(isspace(c) && c != EOF);

    if(c == EOF)
        return (NULL);
    while(!isspace(c) && c != EOF) {
        if(--bufsiz < 1) {
            ungetc(c, ifp);
            break;
        }
        *bp++ = c;
        c = getcic(ifp);
    }
    if(isspace(c))
        ungetc(c, ifp);
    *bp = '\0';
    return (buff);
}

/* read a floating point number and place its value into x */
int
getdblic(double *x, char *buff, int bufsiz, FILE * ifp)
{
    if(getwordic(buff, bufsiz, ifp) == NULL)
        return (EOF);
    sscanf(buff, "%lf", x);
    return (0);
}

/* read an int and place its value into i */
int
getintic(int *i, char *buff, int bufsiz, FILE * ifp)
{
    if(getwordic(buff, bufsiz, ifp) == NULL)
        return (EOF);
    sscanf(buff, "%d", i);
    return (0);
}

/* read an unsigned int and place its value into u */
int
getunsignedic(unsigned *u, char *buff, int bufsiz, FILE * ifp)
{
    if(getwordic(buff, bufsiz, ifp) == NULL)
        return (EOF);
    sscanf(buff, "%u", u);
    return (0);
}
