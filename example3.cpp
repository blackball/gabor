/**
 * A usage example for dense gabor feature.
 *
 * @blackball <bugway@gmail.com>
 */

#include "gabor.h"

#include <cv.h>
#include <highgui.h>

static void _test() {
  int image_w, image_ws, image_h, step_x, step_y;
  int gabor_len;
  gabor_real *feat_vec;
  struct GaborSetting setting;
  IplImage *image = cvLoadImage("008A12.jpg", 0);
  if (!image) return ;

  image_w = image->width;
  image_ws = image->widthStep;
  image_h = image->height;
  step_x = 5;
  step_y = 5;

  gabor_init(&setting, image_w, image_ws, image_h, step_x, step_y);
  gabor_create(&setting);
  gabor_len = gabor_length(&setting);
  feat_vec = (gabor_real*)calloc(gabor_len, sizeof(gabor_real));

  if (!feat_vec) {
    gabor_destroy(&setting);
    return ;
  }

  gabor_extract(&setting, (unsigned char *)image->imageData, feat_vec, GABOR_MAG);

  /*
    do something with feat_vec
  */

  gabor_destroy(&setting);
  cvReleaseImage(&image);
  free(feat_vec);
  
  return ;
}

int main(int argc, char *argv[])
{
  _test();
  getchar();
  return 0;
}
