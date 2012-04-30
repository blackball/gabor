/**
 * A clean and efficient version of Dense gabor feature 
 * based on OpenCV.
 *
 * @TODO will make it lib-free soon.
 *
 * @blackball @date 
 */

#ifndef GABOR_2_H
#define GABOR_2_H

#ifdef __cplusplus
#define EXTERN_BEGIN extern "C" { 
#define EXTERN_END   }
#else
#define EXTERN_BEGIN
#define EXTERN_END
#endif

/* The image size should be know, and we make kernel */
#define IMAGE_WIDTH  84
#define IMAGE_HEIGHT 96
/**
 * kernel_sz = (IMAGE_WIDTH, IMAGE_HEIGHT)
 */
void gabor_makeKernl(kernel_sz, scale, orientation, kernel_vec) {
  make_kernel
}


typedef double gabor_real;
struct GaborSetting {
  int kernel_size;
  double psi;
  double delta;
  double lambda;
  double sigma;
  double x_theta, y_theta;
};

struct GaborFilterBank {
  
};

void gabor_createKernel();
int  gabor_length();
void gabor_extract(const unsigned char* image_data,
                   struct Gabor *...,
                   gabor_real *feat_vec);

#endif 
