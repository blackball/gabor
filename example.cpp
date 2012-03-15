/* unit test */

#include "gabor.h"

static void test() {

  /* image should be gray-scale */
  IMAGE_TYPE *img = getImage();
  int
      imgW = 640,
      imgH = 480,
      gabor_stepW = 10,
      gabor_stepH = 10,
      gabor_len = 0;
  
  float *gabor_vec = NULL;
  
  /* init gabor and return the feature vector length */
  gabor_len = gabor_feature_len(img->width, img->height,
                                gabor_stepW, gabor_stepH);
  gabor_vec = new float[gabor_len];

  /* get feature vector */
  gabor_vec = extract_gabor_features(img, gabor_vec);

  delete gabor_vec;
}


int main(int argc, char *argv[])
{
  test();
  return 0;
}
