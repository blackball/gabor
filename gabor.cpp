
#include "mygabor.h"
#include "extract_gabor.h"

int _imgW = 0;
int _imgH = 0;
int _oritentation = 8;
int _scales = 4;
int _stepH = 0,  _stepW = 0;
int init_flag = 0;

static void _init(int img_w, int img_h, int stepW, int stepH)
{
  _imgW = img_w;
  _imgH = img_h;
  _stepW = stepW;
  _stepH = stepH;
  init_flag = 1;
}

int gabor_feature_len(int img_w, int img_h, int stepW, int stepH)
{
  _init(img_w, img_h, stepW, stepH);
  int extraW = 0, extraH = 0;
  if(img_w % stepW == 0)
  {
    extraW -= 1;
  }
  if(img_h % stepH == 0)
  {
    extraH -= 1;
  }
  return _oritentation * _scales * (img_w/stepW + extraW) * (img_h/stepW + extraH);
}

float* extract_gabor_features(IplImage *gray_img, float *vec)
{
  if(init_flag == 0)
  {
    fprintf(stderr, "Error: Gabor is not initialized!\n");
    exit(-1);
  }
  IplImage *reImg = cvCreateImage(cvGetSize(gray_img), IPL_DEPTH_32F, 1);
  unsigned char *gray_img_data = (unsigned char*)gray_img->imageData;
  float *reimg_data = (float*)reImg->imageData;
  int vec_idx = 0;
  for (int i = 0; i < _scales; ++i)
  {
    for (int j = 0; j < _oritentation; ++j)
    {
      CvGabor gabor(j,i);
      gabor.conv_img(gray_img, reImg, CV_GABOR_MAG);
      // extract by step
      for(int y = _stepH; y < _imgH; y += _stepH)
      {
        for(int x = _stepW; x < _imgW; x += _stepW)
        {
          vec[vec_idx] = reimg_data[y*_imgW + x];
          ++ vec_idx;
        }
      }
			
    }
  }

  cvReleaseImage(&reImg);

  return vec;
}
