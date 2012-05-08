#include "gabor4/gabor.h"
#include <fftw3.h>
#include <cv.h>
#include <highgui.h>

#include <Windows.h>
/*
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
	

	gabor_destroy(&setting);
	cvReleaseImage(&image);
	return ;
}
*/
/*
static void _test_mag() { 
	IplImage *img = cvLoadImage("huo.jpg", 0);
	struct GaborSetting setting;
	gabor_init(&setting, img->width, img->widthStep, img->height);
	gabor_create(&setting);
	for (int i = 0; i < 32; ++i) 
	gabor_on_iplImage(&setting, img, i);
	cvReleaseImage(&img);
	gabor_destroy(&setting);
}
*/

static void printM(double *m, int rows, int cols) {
	for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < cols; ++c) {
			printf("%f,", m[r*cols + c]);
		}
		printf("\n");
	}
}

static void _test_small() {
	double datas[3*3] = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 9);
	fftw_plan plan_f = fftw_plan_dft_r2c_2d(3, 3, datas, out, FFTW_ESTIMATE);
	fftw_plan plan_b = fftw_plan_dft_c2r_2d(3, 3, out, datas, FFTW_ESTIMATE);

	fftw_execute(plan_f);
	fftw_execute(plan_b);

	printM(datas, 3, 3);
}

static void _test4() {
	IplImage *img = cvLoadImage("huo.jpg", 0);
	double *vec;
	int i,j;
	struct gabor_setting setting;
	gabor_init(&setting);
	gabor_create(&setting, img->width, img->height);
	vec = (double*)malloc(sizeof(double) * img->width * img->height);
	cvNamedWindow("test", 1);

	for (int k = 0; k < 40; ++k) {
		gabor_test(&setting, k, (const unsigned char*)img->imageData, img->width, img->widthStep, img->height, vec);
	
		int tmp = 0;
 		for (i = 0; i < img->height; ++i) {
			for (j = 0; j < img->width; ++j) {
				tmp = cvRound(vec[i*img->width + j] * 255.0);
				img->imageData[i*img->widthStep + j] = tmp;
			}
		}	
	cvShowImage("test", img);
	cvWaitKey(0);
	}

	cvDestroyWindow("test");
	cvReleaseImage(&img);
	gabor_destroy(&setting);
}
int main(int argc, char *argv[])
{
	// _test();
	// _test_mag();
	// _test_small();
	_test4();
	getchar();
	return 0;
}