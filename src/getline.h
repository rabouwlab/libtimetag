/* The getline function is in the public domain.
*  Authors: Will Hartung, Antti Haapala
*/

#ifndef GETLINE_H
#define GETLINE_H

#if defined(_WIN32) || defined(_WIN64)

#include <stdio.h>
#include <stdint.h>

intptr_t getline(char **lineptr, size_t *n, FILE *stream);

#endif

#endif // GETLINE_H
