/*
 * Fast Fourier Transform (implementation)
 *
 * Daniel Shved, MIPT, 2010.
 * danshved [at] gmail.com
 */

#include "dft.h"
#include <math.h>
#include <complex>
#include <algorithm>
#include <stack>
using namespace std;

#define PI 3.14159265

/*
 * "Butterfly" transform.
 */
inline void butterfly(comp &x, comp &y, comp w)
{
   comp p = x, q = y*w;
   x = p + q;
   y = p - q;
}

/*
 * Series of butterfly transforms required by the FFT algorithm.
 */
inline void mass_butterfly(comp *array, int size, comp w)
{
   comp power(1.0, 0.0);
   int n = size/2;
   
   for(int i = 0; i < n; i++) {
      butterfly(array[i], array[i+n], power);
      power *= w;
   }
}

/*
 * Given a number ``x'' returns the number which has the same bits as ``x'',
 * but in the reverse order
 */
inline unsigned int backwards(unsigned int x, int length)
{
   unsigned int result = 0;
   unsigned int bit = 1u;
   unsigned int reverse = 1u<<(length-1);
   for(int i = 0; i < length && x != 0; i++) {
      if(x & bit) {
         result |= reverse;
         x &= ~bit;
      }
      bit <<= 1;
      reverse >>= 1;
   }
   return result;
}

/*
 * Moves elements of the array as required by the iterative FFT implementation.
 * ``size'' must be a power of 2.
 */
static void reposition(comp *array, int size)
{
   // Determine the bit length
   int length = 0;
   while(1u << length < (unsigned int)size)
      length++;

   // Swap elements at positions k and reverse(k)
   for(int i = 0; i < size; i++) {
      int j = backwards(i, length);
      if(i <= j)
         swap(array[i], array[j]);
   }
}

/*
 * Does the Discrete Fourier Transform.  Takes time O(size * log(size)).
 * ``size'' must be a power of 2.
 */
void fourier_transform(comp *array, int size)
{
   // Arrange numbers in a convenient order
   reposition(array, size);

   // Prepare roots of unity for every step
   int step;
   comp root = exp(comp(0.0, 2.0*PI/size));
   stack<comp> roots;
   for(step=size; step != 1; step /= 2) {
      roots.push(root);
      root *= root;
   }

   // Do lots of butterfly transforms
   for(step = 2; step <= size; step *= 2) {
      root = roots.top();
      roots.pop();
      for(int i = 0; i < size; i += step)   
         mass_butterfly(array + i, step, root);
   }
}

/*
 * The inverse DFT.
 */
void inverse_fourier_transform(comp *array, int size)
{
   conjugate(array, size);
   fourier_transform(array, size);
   conjugate(array, size);
   for(int i = 0; i < size; i++)
      array[i] = array[i] / (double)size;
}

/*
 * Replaces every element of the vector by its complex conjugate.
 */
void conjugate(comp *array, int size)
{
   for(int i = 0; i < size; i++)
      array[i] = conj(array[i]);
}

/*
 * Multiplies two vectors element by element.
 */
void multiply(comp *arr1, comp *arr2, comp *result, int size)
{
   for(int i = 0; i < size; i++)
      result[i] = arr1[i] * arr2[i];
}

/*
 * Finds the convolution of two vectors (the product of two polynomials, given
 * that the result has power less than ``size'').  ``size'' must be a power of
 * 2.
 */
void convolution(comp *arr1, comp *arr2, comp *result, int size)
{
   fourier_transform(arr1, size);
   fourier_transform(arr2, size);
   multiply(arr1, arr2, result, size);
   inverse_fourier_transform(result, size);
}