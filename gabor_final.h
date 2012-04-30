/**
 * Last version, I misunderstood the dft process when used in convloution.
 * 
 * Here is the final correct version, and for the sake of effcientcy,
 * I made this version configurable.
 */

#ifndef GABOR_H
#define GABOR_H

#ifdef __cplusplus
#define EXTERN_BEGIN extern "C" { 
#define EXTERN_END   }
#else
#define EXTERN_BEGIN
#define EXTERN_END
#endif

EXTERN_BEGIN

typedef double gabor_real;

/**
 * The maximum orientations and scales.
 */
#define GABOR_ORIENTATION_NUM_MAX 8
#define GABOR_SCALE_NUM_MAX 4

/**
 * You need to define a macro to select the mode
 * GABOR_REAL -- just real part convolusion
 * GABOR_IMAG -- just imaginary part convolution
 * GABOR_MAG  -- magnitude ...
 */
#define GABOR_REAL


/**
 * Give all gabor configuration. 
 */
struct GaborBank;
struct GaborSetting {
  int orientation_num;
  int scale_num;

  int image_w;
  int image_ws; /* 'cause I use DFT here, the kernel */
  int image_h; /* size is the same with image size */
  
  int step_x; /* feature sampling step in x-direction */
  int step_y; /* feature sampling step in x-direction */

  int orientations[GABOR_ORIENTATION_NUM_MAX]; /* <num> * PI/8, so num could be [0,7]*/
  int scales[GABOR_SCALE_NUM_MAX]; /* <num> could be [0, 4] */
  /* filter bank */
  struct GaborBank* bank;
};


/**
 * First you need initialize the gabor setting.
 * the function below  is an example, write yours
 * when you're in a different situations.
 *
 * @setting empty gabor setting
 */
void gabor_init(struct GaborSetting *setting,
                int image_w, int image_ws, int image_h,
                int step_x, int step_y);

/**
 * create gabor bank using given setting.
 *
 * @setting filled gabor setting
 */
void gabor_create(struct GaborSetting *setting);

/**
 * Get the gabor feature vector length.
 * 
 * @setting Gabor setting, should be filled before.
 */
int gabor_length(const struct GaborSetting *setting);

/**
 * Extract dense gabor feature.
 *
 */
void gabor_extract(const struct GaborSetting *setting,
                   const unsigned char *image_data,
                   gabor_real *feat_vec, 
                   int gabor_type);

/**
 * Destroy gabor bank.
 */
void gabor_destroy(struct GaborSetting *setting);

EXTERN_END


#endif
