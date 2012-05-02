/**
 * Dense gabor feature, effective and customizable.
 *
 * @blackball @date 5/2/2012
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

/* max number of orientations and scales */
#define GABOR_ORIENTATION_NUM_MAX 8
#define GABOR_SCALE_NUM_MAX 4

/**
 * Before you do the job, 
 * you need to define the Gabor feature
 * type by using #define GABOR_XXX, the MAG
 * one was set to be the default.
 */
/* #define GABOR_REAL */
/* #define GABOR_IMAG */
#define GABOR_MAG

struct GaborBank;
struct GaborSetting {
  int image_w;  /* width */
  int image_ws; /* width step or stride */
  int image_h;  /* height */

  int block_w;  /* block window width */
  int bloch_h;  /* block window height */
  int step_x;   /* step size at x-direction */
  int step_y;   /* step size at x-direction */
  
  int orientation_num;
  int scale_num;
  int orientations[ GABOR_ORIENTATION_NUM_MAX ];
  int scales[ GABOR_SCALE_NUM_MAX ];
  
  struct GaborBank *bank; /* filter banks */
};

/**
 * Initialize GaborSetting, you need define your
 * own setting function, or modify this one.
 *
 * @setting not NULL GaborSetting
 */
void gabor_init(struct GaborSetting *setting);

/**
 * Create Gabor filter bank.
 *
  @setting initialized GaborSetting 
 */
void gabor_create(struct GaborSetting *setting);

/**
 * Get the dense feature length.
 * @setting initialized GaborSetting
 */
int gabor_length(const struct GaborSetting *setting);

/**
 * Extract dense gabor feature from gray-scale
 * image_data into feat_vec;
 *
 */
void gabor_extract(const struct GaborSetting *setting,
                   const unsigned char *image_data,
                   gabor_real *feat*vec);

/**
 * release the gabor bank
 */
void gabor_destroy(struct GaborSetting *setting);

EXTERN_END

#ifdef TEST
/* A usgae example */
{
  struct GaborSetting setting;
  gabor_setting(&setting);
  gabor_create(&setting);
  int feat_len = gabor_length(&setting);
  gabor_real *feat_vec = (gabor_real*)malloc(sizeof(gabor_real) * feat_len);
  ImageType *img = getNextImage(...);
  gabor_extract(&setting, img->image_data, gabor_real);
  /* ... */
  gabor_destroy(&setting);
}
#endif 

#endif
