/**
 * Dense gabor feature based on fftw3.
 *
 * @author blackball (bugway@gmail.com)
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

#define GABOR_SCALE_NUM_MAX 5
#define GABOR_ORIENTATION_NUM_MAX 8

EXTERN_BEGIN

struct gabor_setting {
    int scale_num;
    int orientation_num;
    int scales[ GABOR_SCALE_NUM_MAX ];
    int orientations[GABOR_ORIENTATION_NUM_MAX ];
    
    /* gabor filter bank */
    void *bank;
};

/*** basic gabor api to manage filter bank ****/
int gabor_init(struct gabor_setting *setting);
int gabor_create(struct gabor_setting *setting, int img_w, int img_h);
void gabor_destroy(struct gabor_setting *setting);

/*** Dense gabor feature API ***/
int gabor_length(int image_w,
                 int image_h,
                 int block_w,
                 int block_h,
                 int step_x,
                 int step_y);

int gabor_extract(const struct gabor_setting *setting,
                  const unsigned char *image_data,
                  int image_w,
                  int image_ws,
                  int image_h,
                  double *feat_vec,
                  int feat_len);

void gabor_test(struct gabor_setting *setting, int fidx,
	const unsigned char *image_data, int w, int ws, int h, double *vec);

EXTERN_END

#endif
