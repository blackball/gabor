/**
 * Dense gabor feature implemented using FFTW3.
 *
 * @blackball<bugway@gmail.com>
 */

#ifndef GABOR_FFTW_H
#define GABOR_FFTW_H

#ifdef __cplusplus
#define EXTERN_BEGIN extern "C" { 
#define EXTERN_END   }
#else
#define EXTERN_BEGIN
#define EXTERN_END
#endif

EXTERN_BEGIN

typedef double gabor_real;
#define GABOR_ORIENTATION_NUM_MAX 8
#define GABOR_SCALE_NUM_MAX 4
/**
 * Give all gabor configuration.
 * 
 */
struct GaborBank;
struct GaborSetting {
  int orientation_num;
  int scale_num;

  int kernel_w;
  int kernel_h;

  int image_w; /* 'cause I use DFT here, the kernel */
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
void gabor_init(struct GaborSetting *setting);

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
                   gabor_real *feat_vec);

/**
 * Destroy gabor bank.
 */
void gabor_destroy(struct GaborSetting *setting);

EXTERN_END

#if 0
/* A usage example */
{
  struct GaborSetting setting;
  ImageType *frame = getNextFrame();

  fill_setting( &setting, frame->width, frame->height );
  
  int gabor_len = gabor_length( &setting );
  gabor_real *feat_vec = ( gabor_real* )malloc( sizeof(gabor_real) * gabor_len );
  /* create gabor bank */
  gabor_init(&setting);
  
  while( frame = getNextFrame() ) {
    gabor_extract( &setting, frame->imageData, feat_vec );
    /* do something with feat_vec */
  }

  gabor_destroy( &setting );
  free( feat_vec );  
}
#endif

#endif
