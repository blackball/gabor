/**
 * Dense gabor feature implemented using FFTW3.
 * REAME file has explained why I did in this way. 
 *
 * @blackball(bugway@gmail.com)
 */

#include "gabor3.h"
#include "fftw3.h"
/*****/
#include <stdlib.h> /* malloc */
#include <string.h> /* memset */
#include <math.h>

#define GABOR_PI 3.1415926535897932384626433832795

EXTERN_BEGIN

struct GaborBank {
  int optimal_sz; /* optimal size for fft buffer */

  int filter_num;
  int filter_len;

  void *pMem;
  fftw_complex **bank_dfts;
  fftw_complex *complex_buffer;
  double *real_buffer;
  
  /* using new-array interface */
  fftw_plan plan_forward;
  fftw_plan plan_backward;
};

/**
 * Get optimal FFT size
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
  setting->block_w = 16;
  setting->block_h = 16;
  setting->step_x = step_x;
  setting->step_y = step_y;

  for (i = 0; i < 8; ++i)
    setting->orientations[i] = i;
  for (i = 0; i < 4; ++i)
    setting->scales[i] = i - 1;

  return ;
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
static void _gabor_mk_kernel(int iMu, int iNu,
                             double *real,
                             double *imag,
                             int width, int height) {
  int x,y,i,j;
  
  double Kmax = GABOR_PI / 2;
  double dFreq = sqrt(2.0);
  double K = Kmax / pow(dFreq, (double)iNu);
  double sigma = 2 * GABOR_PI;
  double phi = iMu * (GABOR_PI/8);
  double dTemp1, dTemp2, dTemp3;

  double
      K2 = pow(K, 2),
      sigma2 = pow(sigma, 2),
      k2divsigma2 = K2 / sigma2,
      k2div2sigma2 = K2 / (2*sigma2),
      kcosphi = K * cos(phi),
      ksinphi = K * sin(phi),
      exphalfsigma2 = exp(-sigma2 / 2);

  /* for L1 normalization */
  long double rsum = .0, isum = .0;

  for (i = 0; i < height; ++i) {
    y = i - (height - 1) / 2;
    for (j = 0; j < width; ++j) {
      x = j - (width - 1) / 2;
      dTemp1 = k2divsigma2 * exp( -(x*x+y*y) * k2div2sigma2 );
      dTemp3 = kcosphi * x + ksinphi * y;
      dTemp2 = cos( dTemp3 ) - exphalfsigma2;
      dTemp3 = sin( dTemp3 );

      dTemp2 = dTemp1 * dTemp2;
      rsum += dTemp2;
      real[i*width + j] = dTemp2; 
      dTemp2 = dTemp1 * dTemp3;
      imag[i*width + j] = dTemp2;
      isum += dTemp2;
    }
  }

  /* normalize */
  if (rsum != 0) {
    rsum = 1.0 / rsum;
    for (i = 0; i < height; ++i) {
      for (j = 0; j < width; ++j)
        real[i*width + j] *= rsum;
    }
  }

  if (isum != 0) {
    isum = 1.0 / isum;
    for (i = 0; i < height; ++i) {
      for (j = 0; j < width; ++j)
        imag[i*width + j] *= isum; 
    }
  }
}

/**
 * Shift the kernel center to the corners.
 *
 * @note no need to shift (and pad) the center when in a classic convolution
 *
 * @kernel generated kernel, maybe smaller then the image size
 * @filter image size FFT filter, should be a square in case of a 'unexpected shift'
 */
static void _gabor_pad_kernel(double *kernel, int kw, int kh,
                              double *filter, int fw, int fh){
  

  int ki,kj,fi,fj;
  int center_x = (kw-1) / 2, center_y = (kh-1) / 2;
  double *pKernel = kernel, *pFilter = filter;

  /* kernel size should not bigger than buffer size */
  if (kw > fw || kh > fh) return ;

  for (fi = 0, ki = center_y; ki < kh; ++ki, ++fi) {
    for (fj = 0, kj = center_x; kj < kw; ++kj, ++fj) {
      pFilter[fi*fw + fj] = pKernel[ki*kw + kj];
    }
  }

  for (fi = fh-1, ki = center_y - 1; ki >= 0; --ki, --fi) {
    for (fj = 0, kj = center_x + 1; kj < kw; ++kj, ++fj) {
      pFilter[fi*fw + fj] = pKernel[ki*kw + kj];
    }
  }

  for (fi = 0, ki = center_y + 1; ki < kh; ++ki, ++fi) {
    for (fj = fw - 1, kj = center_x - 1; kj >= 0; --kj, --fj) {
      pFilter[fi*fw + fj] = pKernel[ki*kw + kj];
    }
  }

  for (fi = fh - 1, ki = center_y -1; ki >= 0; --ki, --fi) {
    for (fj = fw - 1, kj = center_x -1 ; kj >= 0; --kj, --fj) {
      pFilter[fi*fw + fj] = pKernel[ki*kw + kj];
    }
  }
}


/**
 * create gabor filter bank, and their DFT represent.
 */
void gabor_create(struct GaborSetting *setting) {

  int i,j,k,sz, opt_sz;
  struct GaborBank *pBank;
  fftw_complex *tmp;
  double *imag_buffer;
  
  setting->bank = (struct GaborBank*)malloc(sizeof(struct GaborBank));
  pBank = setting->bank;

  pBank->filter_num = setting->scale_num * setting->orientation_num;

  i = _gabor_getoptimalsz(setting->image_w);
  j = _gabor_getoptimalsz(setting->image_h);

  opt_sz = i > j ? i : j;
  pBank->optimal_sz = opt_sz;

  pBank->filter_len = opt_sz * opt_sz;
  
  /* allocate for filters and plans */
  sz  = sizeof(fftw_complex*) * pBank->filter_num * 2;
  sz += sizeof(fftw_complex)  * (pBank->filter_num) * pBank->filter_len * 2;
  sz += sizeof(fftw_complex) * pBank->filter_len; /* for complex_buffer */
  sz += sizeof(double) * pBank->filter_len; /* for real_buffer */
  
  pBank->pMem = fftw_malloc(sz);
  
  if (!pBank->pMem) return ;

  memset(pBank->pMem, 0, sz);
  
  pBank->bank_dfts = (fftw_complex**)(pBank->pMem);
  tmp = (fftw_complex*)(pBank->bank_dfts + pBank->filter_num * 2);

  for (i = 0; i < pBank->filter_num * 2; i += 2) {
    *(pBank->bank_dfts+i) = tmp;
    tmp += pBank->filter_len; /* real part */
    *(pBank->bank_dfts + i + 1) = tmp;
    tmp += pBank->filter_len; /* imag part */
  }
  
  pBank->idft_buffer = tmp; /* complex_buffer */
  tmp += pBank->filter_len;
  pBank->real_buffer = (double*)tmp; /* real_buffer */

  /* create plans r2c and c2r */
  pBank->plan_forward = fftw_plan_dft_r2c_2d(opt_sz,
                                             opt_sz,
                                             pBank->real_buffer,
                                             pBank->complex_buffer,
                                             FFTW_MEASURE);
  pBank->plan_backward = fftw_plan_dft_c2r_2d(opt_sz,
                                              opt_sz,
                                              pBank->complex_buffer,
                                              pBank->real_buffer,
                                              FFTW_MEASURE);
  
  imag_buffer = (double*)calloc(pBank->filter_len, sizeof(double));
  pad_buffer =  (double*)calloc(pBank->filter_len, sizeof(double));

  if (!imag_buffer || !pad_buffer) return ;

  k = 0;
  for (i = 0; i < setting->scale_num; ++i)
    for (j = 0; j < setting->orientation_num; ++j) {

      _gabor_mk_kernel(setting->orientations[j],
                       setting->scales[i],
                       pBank->real_buffer, /* real part kernel */
                       pBank->imag_buffer, /* imag part kernel */
                       opt_sz,
                       opt_sz);
      /* shift the kenerl to the 0-padding image size filter */
      _gabor_pad_kernel(pBank->real_buffer, opt_sz, opt_sz,
                        pad_buffer, opt_sz, opt_sz);
      /* real part */
      fftw_execute_dft_r2c(pBank->plan_forward,
                           pBank->pad_buffer,
                           pBank->bank_dfts[k]);
      /* shift the kenerl to the 0-padding image size filter */
      _gabor_pad_kernel(pBank->imag_buffer, opt_sz, opt_sz,
                        pad_buffer, opt_sz, opt_sz); 
      /* imag part */
      fftw_execute_dft_r2c(pBank->plan_forward,
                           pBank->pad_buffer,
                           pBank->bank_dfts[k+1]);
      k += 2;
    }
  
  free(imag_buffer);
  free(pad_buffer);
  
  return ;
}

int gabor_length(const struct GaborSetting *setting) {
  /* considering the block representation */
  return 0; 
}

static void _gabor_norm_real(fftw_complex *t, int len, gabor_real scale) {
  int i;
  fftw_complex *p = t;
  gabor_real rmid, rmax, rmin;

  rmax = rmin = p[0][0];
  for (i = len - 1; i--; p++) {
    rmin = rmin > p[0][0] ? p[0][0] : rmin;
    rmax = rmax < p[0][0] ? p[0][0] : rmax;
  }

  rmid = scale / (rmax - rmin);
  p = t;
  
  if (rmid != 0)
    for (i = len - 1; i--; p++) {
      p[0][0] -= rmin;
      p[0][0] *= rmid;
    }
}

static void _gabor_norm_imag(fftw_complex *t, int len, gabor_real scale) {
  int i;
  fftw_complex *p = t;
  gabor_real imid, imax, imin;

  imax = imin = p[0][1];
  for (i = len - 1; i--; p++) {
    imin = imin > p[0][1] ? p[0][1] : imin;
    imax = imax < p[0][1] ? p[0][1] : imax;
  }

  imid = scale / (imax - imin);
  p = t;
  
  if (imid != 0)
    for (i = len - 1; i--; p++) {
      p[0][1] -= imin;
      p[0][1] *= imid;
    }
}

static void _gabor_norm_mag(fftw_complex *t, int len, gabor_real scale) {
  int i;
  fftw_complex *p = t;
  gabor_real rmid, imid, rmax, rmin, imax, imin;

  rmax = rmin = p[0][0];
  imax = imin = p[0][1];
  for (i = len - 1; i--; p++) {
    rmin = rmin > p[0][0] ? p[0][0] : rmin;
    rmax = rmax < p[0][0] ? p[0][0] : rmax;
    imin = imin > p[0][1] ? p[0][1] : imin;
    imax = imax < p[0][1] ? p[0][1] : imax;
  }

  rmid = scale / (rmax - rmin);
  imid = scale / (imax - imin);

  p = t;
  if (rmid != 0)
    for (i = len - 1; i--; p++) {
      p[0][0] -= rmin;
      p[0][0] *= rmid;
    }
  
  p = t;
  if (imid != 0)
    for (i = len - 1; i--; p++) {
      p[0][1] -= imin;
      p[0][1] *= imid;
    }
}

#if 0
#include <cv.h>
#include <highgui.h>
static void _gabor_test(fftw_complex *t, int opt_sz, int imgw, int imgh) {
  int i,j;
  IplImage *img = cvCreateImage(cvSize(imgw, imgh), 8, 1);
	
  for (i = 0; i < imgh; ++i) {
    for (j = 0; j < imgw; ++j) {
      img->imageData[i*img->widthStep + j] = t[i*opt_sz + j][0];//sqrt(t[i*opt_sz + j][0] * t[i*opt_sz + j][0] + t[i*opt_sz + j][1] * t[i*opt_sz + j][1]);
    }
  }

  cvNamedWindow("test", 1);
  cvShowImage("test", img);
  cvWaitKey(0);
  cvDestroyWindow("test");
  cvReleaseImage( &img );
}
#endif

void gabor_extract(const struct GaborSetting *setting,
                   const unsigned char *image_data,
                   gabor_real *feat_vec,
                   int gabor_type) {
  /* extract dense gabor feature */
  struct GaborBank *pBank = setting->bank;
  int filter_num = pBank->filter_num;
  int image_w = setting->image_w;
  int image_ws = setting->image_ws;
  int image_h = setting->image_h;
  int step_x = setting->step_x;
  int step_y = setting->step_y;
  int opt_sz = setting->bank->opt_sz;
  
  int i, j, _x, _y;
  gabor_real *pFeat;

  fftw_complex
      *in_buffer = pBank->in_buffer,
      *dft_buffer = pBank->dft_buffer,
      *idft_buffer = pBank->idft_buffer,
      *pCurrKernel,
      tmp0, tmp1;
  
  /* fill in image data */
  
  int filter_len = pBank->filter_len;
  memset(in_buffer, 0, sizeof(fftw_complex) * filter_len);

  for (i = 0; i < image_h; ++i) {
    for (j = 0; j < image_w; ++j) {
      in_buffer[i*opt_sz + j][0] = image_data[i * image_ws + j];
      in_buffer[i*opt_sz + j][1] = 0.0;
    }
  }
  
  /* DFT on image data and store image's dft data in dft_buffer */
  fftw_execute_dft(pBank->plan_forward,
                   in_buffer,
                   dft_buffer);
  
  pFeat = feat_vec;
  
  for (i = 0; i < filter_num; ++i)  {

    /* multiply the dft vectors */
    pCurrKernel = pBank->bank_dfts[i];
    for (j = 0; j < filter_len; ++j) {
      tmp0[0] = dft_buffer [j][0];
      tmp0[1] = dft_buffer [j][1];
      tmp1[0] = pCurrKernel[j][0];
      tmp1[1] = pCurrKernel[j][1];
      in_buffer[j][0] = (tmp0[0] * tmp1[0] - tmp0[1] * tmp1[1]) ;
      in_buffer[j][1] = (tmp0[0] * tmp1[1] + tmp0[1] * tmp1[0]) ;
    }
    
    /* inverse dft */
    fftw_execute_dft(pBank->plan_backward,
                     in_buffer,
                     idft_buffer);
    
    /* extract MAG features */
    switch(gabor_type) {
      
      case GABOR_REAL:
        _gabor_norm_real(idft_buffer, filter_len, 255.0);

        for (_y = 0; _y < image_h; _y += step_y) {
          for (_x = 0; _x < image_w; _x += step_x) {
            *pFeat++ = idft_buffer[_y * opt_sz + _x][0];
          }
        }
        break;
        
      case GABOR_IMAG:
        _gabor_norm_imag(idft_buffer, filter_len, 255.0);
        

        for (_y = 0; _y < image_h; _y += step_y) {
          for (_x = 0; _x < image_w; _x += step_x) {
            *pFeat++ = idft_buffer[_y * opt_sz + _x][1];
          }
        }
        break;
        
      case GABOR_MAG:
        _gabor_norm_mag(idft_buffer, filter_len, 255.0);
        
#if 0
        /***TEST***/
        _gabor_test(idft_buffer, pBank->opt_sz, setting->image_w, setting->image_h);
#endif 
        for (_y = 0; _y < image_h; _y += step_y) {
          for (_x = 0; _x < image_w; _x += step_x) {
            tmp0[0] = idft_buffer[_y * opt_sz + _x][0];
            tmp0[1] = idft_buffer[_y * opt_sz + _x][1];
            /* @TODO you maybe need to do a min-max normalization on pFeat */
            *pFeat++ = sqrt(tmp0[0]*tmp0[0] + tmp0[1]*tmp0[1]);
          }
        }
        break;
        
      default:
        _gabor_norm_real(idft_buffer, filter_len, 1.0/255);
        for (_y = 0; _y < image_h; _y += step_y) {
          for (_x = 0; _x < image_w; _x += step_x) {
            *pFeat++ = idft_buffer[_y * opt_sz + _x][0];
          }
        }
    }
  }
  return ;
}

void gabor_destroy(struct GaborSetting *setting) {
  /* free resources in gabor bank */
  if ( setting && setting->bank ) {
    fftw_free( setting->bank->pMem );
    fftw_destroy_plan( setting->bank->plan_forward );
    fftw_destroy_plan( setting->bank->plan_backward );
    free( setting->bank );
    setting->bank = NULL;
  }
  return ;
}

EXTERN_END
