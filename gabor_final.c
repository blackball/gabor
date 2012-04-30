/**
 * Dense gabor feature implemented using FFTW3.
 * REAME file has explained why I did in this way. 
 *
 * @blackball(bugway@gmail.com)
 */

#include "gabor_final.h"
#include "fftw3.h"

#include <stdlib.h>
#include <stding.h>

EXTERN_BEGIN

#define GABOR_PI 3.1415926535897932384626433832795


struct GaborBank {
  int opt_sz; /* optimal size for fft buffer */

  int filter_num;
  int filter_len;

  void *pMem;
  fftw_complex **bank_dfts;
  fftw_complex *in_buffer;
  fftw_complex *dft_buffer;
  fftw_complex *idft_buffer; /* hold inverse*/
  /* using new-array interface */
  fftw_plan plan_forward;
  fftw_plan plan_backward;
};

/**
 * Get optimal DFT size
 */
static int _gabor_getoptimalsz(int n) {
  int newsz = n,m;
  for(;;) {
    m = newsz;
    while((m%2) == 0) m /= 2;
    while((m%3) == 0) m /= 3;
    while((m%5) == 0) m /= 5;
    if (m <= 1) break;
    newsz++;
  }
  return newsz;
}


/**
 * Custom initialize the gabor setting.
 *
 * @notice you may define your own gabor_init for your case.
 *
 * @setting un-initialized gaboe setting
 * @img_w width of image data
 * @img_ws width step of image data (or 'stride', in another word)
 * @img_h height of image data
 * @step_x step size in x-direction
 * @step_y step size in y-direction
 */
void gabor_init(struct GaborSetting *setting,
                int img_w, int img_ws, int img_h,
                int step_x, int step_y) {
  int i;

  setting->orientation_num = 8;
  setting->scale_num = 4;
  setting->image_w = img_w;
  setting->image_ws = img_ws;
  setting->image_h = img_h;
  setting->step_x = step_x;
  setting->step_y = step_y;
  
  for (i = 0; i < 8; ++i)
    setting->orientations[i] = i;
  for (i = 0; i < 4; ++i)
    setting->scales[i] = i - 1;
}

/**
 * Make gabor kernel and transform & align 
 * to be proper to do DFT.
 *
 * @iMu the orientation iMu*PI/8, [0,7]
 * @iNu the scale, [0,4]
 * @vec kernel includeing real and imaginary parts
 * @width kernel width
 * @height kernel height 
 */
static void _gabor_mk_kernel(int iMu,
                             int iNu,
                             fftw_complex *vec,
                             int width,
                             int height) {
  int x,y,i,j;
  
  double Kmax = GABOR_PI / 2;
  double dFreq = sqrt(2.0);
  double K = Kmax / pow(dFreq, (double)iNu);
  double sigma = 2 * GABOR_PI;
  double phi = iMu * (GABOR_PI/8);
  double dTemp0;
  double
      dTemp1,
#if   defined( GABOR_REAL )
      dTemp2, rsum
#elif defined( GABOR_IMAG )
      dTemp3, isum,
#elif defined( GABOR_MAG )
      dTemp2, dTemp3, rsum, isum
#endif
      ;
      
  double
      K2 = pow(K, 2),
      sigma2 = pow(sigma, 2),
      k2divsigma2 = K2 / sigma2,
      k2div2sigma2 = K2 / (2*sigma2),
      kcosphi = K * cos(phi),
      ksinphi = K * sin(phi),
      exphalfsigma2 = exp(-sigma2 / 2);

  rsum = .0; isum = .0;
  for (i = 0; i < height; ++i) {
    y = i - (height - 1) / 2;
    for (j = 0; j < width; ++j) {
      x = j - (width - 1) / 2;
      dTemp1 = k2divsigma2 * exp( -(x*x+y*y) * k2div2sigma2 );
      dTemp0 = kcosphi * x + ksinphi * y;

#if   defined( GABOR_REAL )
      dTemp2 = cos( dTemp0 ) - exphalfsigma2;
      dTemp2 = dTemp1 * dTemp2;
      rsum += dTemp2;
      (vec[i*width + j])[0] = dTemp2;
#elif defined( GABOR_IMAG)
      dTemp3 = sin( dTemp0 );
      dTemp3 = dTemp1 * dTemp3;
      isum += dTemp3;
      ()
#elif defined( GABOR_MAG )
         @TODO
#endif
      dTemp2 = cos( dTemp0 ) - exphalfsigma2;
      dTemp3 = sin( dTemp0 );

      dTemp2 = dTemp1 * dTemp2;
      rsum += dTemp2;
      (vec[i*width + j])[0] = dTemp2; 
      dTemp2 = dTemp1 * dTemp3;
      (vec[i*width + j])[1] = dTemp2;
      isum += dTemp2;
    }
  }

  /* normalize */
  if (rsum != 0) {
    rsum = 1.0 / rsum;
    for (i = 0; i < height; ++i) {
      for (j = 0; j < width; ++j)
        (vec[i*width + j])[0] *= rsum; 
    }
  }

  if (isum != 0) {
    isum = 1.0 / isum;
    for (i = 0; i < height; ++i) {
      for (j = 0; j < width; ++j)
        (vec[i*width + j])[1] *= isum; 
    }
  }
}


EXTERN_END
