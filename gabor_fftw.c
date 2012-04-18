/**
 * Dense gabor feature implemented using FFTW3.
 * REAME file has explained why I did in this way. 
 *
 * @blackball(bugway@gmail.com)
 */

#include "gabor_fftw.h"
#include "fftw3.h"
#include <math.h>

#define GABOR_PI 3.1415926535897932384626433832795

struct GaborBank {
  int filter_num;
  int filter_len;

  void *pMem;
  fftw_complex **bank_dfts;
  fftw_complex *in_buffer;
  fftw_complex *dft_buffer;
  fftw_complex *idft_buffer; /* hold inverse*/
  /* using new-array interface */
  fftw_plan plan_foreward;
  fftw_plan plan_backward;
};

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
  int x,y;
  double dTemp1, dTemp2, dTemp3;
  double Kmax = GABOR_PI / 2;
  double dFreq = sqrt(2.0);
  double K = Kmax / pow(dFreq, (double)iNu);
  double sigma = 2 * GABOR_PI;
  double phi = iMu * (GABOR_PI/8);
  
  double
      K2 = pow(K, 2),
      sigma2 = pow(sigma, 2),
      k2divsigma2 = K2 / sigma2,
      k2div2sigma2 = K2 / (2*sigma2),
      kcosphi = K * cos(phi),
      ksinphi = K * sin(phi),
      exphalfsigma2 = exp(-sigma2 / 2);
  
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      x = i - width / 2;
      y = i - height / 2;
      
      dTemp1 = k2divsigma2 * exp( -(x*x+y*y) * k2div2sigma2 );
      dTemp2 = cos( kcosphi * x + ksinphi * y ) - exphalfsigma2;
      dTemp3 = sin( kcosphi * x + ksinphi * y );

      (vec + i*width + j)->real = dTemp1 * dTemp2;
      (vec + i*width + j)->imag = dTemp1 * dTemp3;
    }
  }
}

/**
 * Initialize gabor filter bank, and their DFT represent.
 */
void gabor_init(struct GaborSetting *setting) {
  struct GaborBan *pBank;
  fftw_complex *tmp;
  
  if (!setting) return ;
  
  setting->bank = (struct GaborBank*)malloc(sizeof(struct GaborBank));
  pBank = setting->bank;
  pBank->filter_num = setting->scale_num * setting->orientation_num;
  pBank->filter_len = setting->kenerl_w * setting->kenerl_h;
  /* allocate for filters and plans */

  int sz = 0;
  sz  = sizeof(fftw_complex*) * pBank->filter_num;
  sz += sizeof(fftw_complex)  * (pBank->filter_num + 3) * pBank->filter_len;
  
  pBank->pMem = fftw_malloc(sz);
  
  if (!pBank->pMem) return ;
  memset(pBank->pMem, 0, sz);
  
  pBank->bank_dfts = (fftw_complex**)pMem;
  tmp = (fftw_complex*)(pBank->bank_dfts + pBank->filter_num);

  for (int i = 0; i < pBank->filter_num; ++i) {
    *(pBank->bank_dfts+i) = tmp;
    tmp += pBank->filter_sz;
  }
  pBank->in_buffer = tmp;
  tmp += pBank->filter_sz;
  pBank->dft_buffer = tmp;
  tmp += pBank->filter_sz;
  pBank->idft_buffer = tmp;

  pBank->plan_foreward = fftw_create_dft_plan_2d(setting->kenerl_w,
                                                 setting->kenerl_h,
                                                 in_buffer,
                                                 dft_buffer,
                                                 FFTW_FOREWARD,
                                                 FFTW_ESTIMATE);
  
  pBank->plan_backward = fftw_create_dft_plan_2d(setting->kenerl_w,
                                                 setting->kenerl_h,
                                                 dft_buffer,
                                                 idft_buffer,
                                                 FFTW_BACKWARD,
                                                 FFTW_ESTIMATE);
  
  for (int i = 0; i < pBank->filter_num; ++i) {
    /* fill bank */
    _gabor_mk_kernel(setting->iMus[i],
                     setting->iNus[i],
                     pBank->in_buffer,
                     setting->kernel_w,
                     setting->kernel_h);
    
    /* get the dft of the kernel */
    fftw_execute_dft(plan_foreward,
                     pBank->in_buffer,
                     pBank->bank_dfts[i]);
  }

  return ;
}

int gabor_length(const struct GaborSetting *setting,
                 int width,
                 int height) {
  return
      (setting->kernel_w / setting->step_x) *
      (setting->kernel_h / setting->step_y) *
      (setting->orientation_num) * (setting->scale_num);
}

void gabor_extract(const struct GaborSetting *setting,
                   const unsigned char *image_data,
                   gabor_real *feat_vec) {
  /* extract dense gabor feature */
  struct GaborBank *pBank = setting->bank;
  int filter_num = pBank->filter_num;
  int kernel_w = setting->kenerl_w;
  int kernel_h = setting->kenerl_h;
  int step_x = setting->step_x;
  int step_y = setting->step_y;
  
  fftw_complex
      *in_buffer = pBank->in_buffer,
      *dft_buffer = pBank->dft_buffer,
      *idft_buffer = pBank->idft_buffer;
  
  /* fill in image data */
  int filter_len = pBank->filter_len;
  for (int i = 0; i < filter_len; ++i) {
    in_buffer[i].real = image_data[i];
    in_buffer[i],imag = .0;
  }
  /* DFT on image data and store image's dft data in dft_buffer */
  fftw_execute_dft(pBank->plan_foreward,
                   in_buffer,
                   dft_buffer);

  gabor_real *pFeat = feat_vec;
  for (int i = 0; i < filter_num; ++i) {
    /* multiply the dft vectors */
    fftw_complex *pCurrKernel = pBank->bank_dfts[i];
    for (int j = 0; j < filter_len; ++j) {
      fftw_complex tmp0 = dft_buffer[j], tmp1 = pCurrKernel[j];
      dft_buffer[j].real = (tmp0.real * tmp1.real - tmp0.imag * tmp1.imag);
      dft_buffer[j].imag = (tmp0.real * tmp1.imag + tmp0.imag * tmp1.real);
    }
    /* inverse dft */
    fftw_execute_dft(pBank->plan_backward,
                     dft_buffer,
                     idft_buffer);
    /* extract MAG features */
    for (int _y = 0; _y < kenerl_h; _y += step_y) {
      for (int _x = 0; _x < kernel_w; _x += step_x) {
        fftw_complex tmp2 = idft_buffer[ _y * kenerl_w + _x ];
        *pFeat++ = sqrt(tmp2.real * tmp2.real + tmp2.imag * tmp2.imag );
      }
    }
  }

  return ;
}

void gabor_destroy(struct GaborSetting *setting) {
  /* free resources in gabor bank */
  if ( setting && setting->bank ) {
    fftw_free( setting->bank->pMem );
    fftw_destroy_plan( setting->bank->plan_foreward );
    fftw_destroy_plan( setting->bank->plan_backword );
    free( setting->bank );
    setting->bank = NULL;
  }
  return ;
}
