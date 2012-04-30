/**
 * Simplified gabor based on OpenCV.
 * @TODO reimplement the OpenCV part for more efficiency
 * @blackball
 */

#include <cv.h>
#include <highgui.h>

/* create kernel bank */

int gabor(const CvMat* src,
          const CvMat* kernel,
          CvMat *dst,
          int type); /*REAL, IMAG, MAG*/

/**
 * Re-implement dense Gabor feathure usig FFTW3. 
 *
 * First, create all gabor kernels' FFT plane, this need
 * to extend the kernel size as the same with the image.
 * Then, for every image, calculate its FFT, multiply the
 * filters' FFT bank, and inverse.
 *
 * Given the image is small(84 * 96), and kernel size(15 * 15) is
 * not small, the FFTW pre-calculation will save some time.
 */


