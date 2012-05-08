/**
 * Dense gabor feature implemented using FFTW3.
 * REAME file has explained why I did in this way. 
 *
 * @blackball(bugway@gmail.com)
 */

#include "gabor.h"
#include "fftw3.h"
/*****/
#include <stdlib.h> /* malloc */
#include <string.h> /* memset */
#include <math.h>

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
 * Get optimal size for FFT
 *
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
 * @setting un-initialized gaboe setting
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
static void _gabor_mk_kernel(int iMu, int iNu, fftw_complex *vec, int width, int height) {
  int x,y,i,j;
  
  double Kmax = GABOR_PI / 2;
  double dFreq = sqrt(2.0);
  double K = Kmax / pow(dFreq, (double)iNu);
  double sigma = 2 * GABOR_PI;
  double phi = iMu * (GABOR_PI/8);
  double dTemp1, dTemp2, dTemp3;

  /**
  * @TODO get optimal kernel size
  */

  double
      K2 = pow(K, 2),
      sigma2 = pow(sigma, 2),
      k2divsigma2 = K2 / sigma2,
      k2div2sigma2 = K2 / (2*sigma2),
      kcosphi = K * cos(phi),
      ksinphi = K * sin(phi),
      exphalfsigma2 = exp(-sigma2 / 2);
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

/**
 * Shift the kernel center to the corners.
 *
 * @TODO Directly shift kernel in generating process.
 */
static void _gabor_pad_kernel(fftw_complex *kernel, int kw, int kh,
                              fftw_complex *filter, int fw, int fh){
  

  int ki,kj,fi,fj;
  int center_x = (kw-1) / 2, center_y = (kh-1) / 2;
  fftw_complex *pKernel = kernel, *pFilter = filter;

  /* kernel size should not bigger than buffer size */
  if (kw > fw || kh > fh) return ;

  for (fi = 0, ki = center_y; ki < kh; ++ki, ++fi) {
    for (fj = 0, kj = center_x; kj < kw; ++kj, ++fj) {
      pFilter[fi*fw + fj][0] = pKernel[ki*kw + kj][0];
      pFilter[fi*fw + fj][1] = pKernel[ki*kw + kj][1];
    }
  }

  for (fi = fh-1, ki = center_y - 1; ki >= 0; --ki, --fi) {
    for (fj = 0, kj = center_x + 1; kj < kw; ++kj, ++fj) {
      pFilter[fi*fw + fj][0] = pKernel[ki*kw + kj][0];
      pFilter[fi*fw + fj][1] = pKernel[ki*kw + kj][1];
    }
  }

  for (fi = 0, ki = center_y + 1; ki < kh; ++ki, ++fi) {
    for (fj = fw - 1, kj = center_x - 1; kj >= 0; --kj, --fj) {
      pFilter[fi*fw + fj][0] = pKernel[ki*kw + kj][0];
      pFilter[fi*fw + fj][1] = pKernel[ki*kw + kj][1];
    }
  }

  for (fi = fh - 1, ki = center_y -1; ki >= 0; --ki, --fi) {
    for (fj = fw - 1, kj = center_x -1 ; kj >= 0; --kj, --fj) {
      pFilter[fi*fw + fj][0] = pKernel[ki*kw + kj][0];
      pFilter[fi*fw + fj][1] = pKernel[ki*kw + kj][1];
    }
  }
}


/**
 * create gabor filter bank, and their DFT represent.
 */
void gabor_create(struct GaborSetting *setting) {

  int i,j,sz;
  struct GaborBank *pBank;
  fftw_complex *tmp;
  
  if (!setting) return ;
  
  setting->bank = (struct GaborBank*)malloc(sizeof(struct GaborBank));
  pBank = setting->bank;

  pBank->filter_num = setting->scale_num * setting->orientation_num;

  i = _gabor_getoptimalsz(setting->image_w);
  j = _gabor_getoptimalsz(setting->image_h);

  i = i > j ? i : j;
  
  pBank->opt_sz = i;

  pBank->filter_len = i*i;
  
  /* allocate for filters and plans */

  sz  = sizeof(fftw_complex*) * pBank->filter_num;
  sz += sizeof(fftw_complex)  * (pBank->filter_num + 3) * pBank->filter_len;
  
  pBank->pMem = fftw_malloc(sz);
  
  if (!pBank->pMem) return ;

  memset(pBank->pMem, 0, sz);
  
  pBank->bank_dfts = (fftw_complex**)(pBank->pMem);
  tmp = (fftw_complex*)(pBank->bank_dfts + pBank->filter_num);

  for (i = 0; i < pBank->filter_num; ++i) {
    *(pBank->bank_dfts+i) = tmp;
    tmp += pBank->filter_len;
  }
  
  pBank->in_buffer = tmp;
  tmp += pBank->filter_len;
  pBank->dft_buffer = tmp;
  tmp += pBank->filter_len;
  pBank->idft_buffer = tmp;

  pBank->plan_forward = fftw_plan_dft_2d(pBank->opt_sz,
                                          pBank->opt_sz,
                                          pBank->in_buffer,
                                          pBank->dft_buffer,
                                          FFTW_FORWARD,
                                          FFTW_ESTIMATE);
  
  pBank->plan_backward = fftw_plan_dft_2d(pBank->opt_sz,
                                          pBank->opt_sz,
                                          pBank->dft_buffer,
                                          pBank->idft_buffer,
                                          FFTW_BACKWARD,
                                          FFTW_ESTIMATE);

  for (i = 0; i < setting->scale_num; ++i)
    for (j = 0; j < setting->orientation_num; ++j) {

      _gabor_mk_kernel(setting->orientations[j],
                       setting->scales[i],
                       pBank->dft_buffer, /* create kernel temporarily stores here */
                       pBank->opt_sz,
                       pBank->opt_sz);

      /* shift the kenerl to the 0-padding image size filter */
      _gabor_pad_kernel(pBank->dft_buffer, pBank->opt_sz, pBank->opt_sz,
                        pBank->in_buffer, pBank->opt_sz, pBank->opt_sz); 

      /* get the dft of the kernel */
	
      fftw_execute_dft(pBank->plan_forward,
                       pBank->in_buffer,
                       pBank->bank_dfts[ i*setting->orientation_num + j ]);
      
    }
  
  return ;
}

int gabor_length(const struct GaborSetting *setting) {

  int opt_w, opt_h, opt_sz, nx, ny;
  opt_w = _gabor_getoptimalsz(setting->image_w);
  opt_h = _gabor_getoptimalsz(setting->image_h);
  
  /* get the bigger one */
  opt_sz = opt_w > opt_h ? opt_w : opt_h;
  
  nx = (opt_sz - opt_sz % setting->step_x) / setting->step_x + 1;
  ny = (opt_sz - opt_sz % setting->step_y) / setting->step_y + 1;
  
  return
      (nx * ny) * (setting->orientation_num) * (setting->scale_num);
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

#if 1
#include <cv.h>
#include <highgui.h>
static void _gabor_test(fftw_complex *t, int opt_sz, int imgw, int imgh) {
	int i,j;
	IplImage *iimg = cvCreateImage(cvSize(imgw, imgh), 8, 1);
	IplImage *rimg = cvCreateImage(cvSize(imgw, imgh), 8, 1);
	for (i = 0; i < imgh; ++i) {
		for (j = 0; j < imgw; ++j) {
			rimg->imageData[i*imgw + j] = cvRound(t[i*opt_sz +j][0]);
			iimg->imageData[i*imgw + j] = cvRound(t[i*opt_sz +j][1]);
		}
	}

	cvNamedWindow("real", 1);
	cvNamedWindow("imag", 1);
	cvShowImage("real", rimg);
	cvShowImage("imag", rimg);
	cvWaitKey(0);
	cvDestroyWindow("real");
	cvDestroyWindow("imag");
	cvReleaseImage( &rimg );
	cvReleaseImage( &iimg );
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
        
#if 1
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
