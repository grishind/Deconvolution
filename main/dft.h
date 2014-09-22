/*
 * Fast Fourier Transform (with complex numbers)
 *
 * Daniel Shved, MIPT, 2010.
 * danshved [at] gmail.com
 */

#ifndef __DFT_H__
#define __DFT_H__

#include <complex>
using namespace std;

typedef complex<double> comp;

// Gets the complex conjugate of every element 
void conjugate(comp *array, int size);

// Multiplies two vectors element by element
void multiply(comp *arr1, comp *arr2, comp *result, int size);

// Finds the convolution of two vectors
void convolution(comp *arr1, comp *arr2, comp *result, int size);

// Discrete fourier transform. size must be a power of 2
void fourier_transform(comp *array, int size);

// Inverse fourier transform. size must be a power of 2
void inverse_fourier_transform(comp *array, int size);

#endif