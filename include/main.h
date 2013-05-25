#ifndef MAIN_H
#define	MAIN_H

#define _USE_MATH_DEFINES

#include "interp.h"
#include "globals.h"
#include "debug.h"
#include "functions.h"
#include "output.h"
#include "border.h"
#include "init.h"
#include "projectile.h"

#include <ctime>
#include <cstring> 
#include <stdexcept>
#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <omp.h>




/* Base64 encoding */
static const std::string base64_chars = 
             "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
             "abcdefghijklmnopqrstuvwxyz"
             "0123456789+/";


static inline bool is_base64(unsigned char c) {
  return (isalnum(c) || (c == '+') || (c == '/'));
}

double truncNdigit(double value, int N);
double truncNdigit(double value, int N);
double Bi(unsigned int n, int i, double x);
double Bj(unsigned int m, int j, double y);
double splineEval(double x, double y, int n, int m, double ** k);
unsigned int factorial(unsigned int n);
std::string base64_encode(unsigned char const* , unsigned int len);


#endif	/* MAIN_H */
