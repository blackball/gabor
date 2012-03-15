

#ifndef _EXTRACT_GABOR_H_
#define _EXTRACT_GABOR_H_

typedef struct _IplImage IplImage;

int gabor_feature_len(int img_w, int img_h, int stepW, int stepH);
float* extract_gabor_features(IplImage *gray_img, float *vec);

#endif
