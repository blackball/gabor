/**
 * implement of dense gabor feature based on FFTW3
 *
 * @author blackball (bugway@gmail.com)
 */

#include "gabor.h"
#include <fftw3.h>
#include <string.h>
#include <assert.h>

#ifndef GABOR_PI
#define GABOR_PI 3.1415926535897932384626433832795
#endif

EXTERN_BEGIN

#include "statics.c"

struct gabor_bank {
    int opt_sz; /* optimal size for filter */

    fftw_plan plan_forward;
    fftw_plan plan_backward;

    void *pmem; /* memory chunk pointer */
    fftw_complex **filters;
    fftw_complex *img_dfted; /* buffer for DFTed image */
    fftw_complex *real_mult; /* multipliaction of real part */
    fftw_complex *imag_mult; /* multiplication of imag part */
    double *buff0;
    double *buff1;
};

int gabor_init(struct gabor_setting *setting) {
    int i;

    setting->scale_num = 5;
    setting->orientation_num = 8;

    for (i = 0; i < setting->scale_num; ++i)
        setting->scales[ i ] = i - 1;
    for (i = 0; i < setting->orientation_num; ++i)
        setting->orientations[ i ] = i;

    setting->bank = 0;

    return 0;
}

void gabor_destroy(struct gabor_setting *setting) {
    if (!setting && !setting->bank) {
        struct gabor_bank *pbank = (struct gabor_bank*)setting->bank;
        fftw_free( pbank->pmem );
        fftw_free( setting->bank );
        setting->bank = NULL;
    }
}

/**
 * create gabor bank filters by initialized setting.
 */
int gabor_create(struct gabor_setting *setting, int img_w, int img_h) {
    int i, j, k;
    int real_num, imag_num, sz, opt_sz, len;
    struct gabor_bank *pbank;
    double *real_kernel_buff, *imag_kernel_buff;
	fftw_complex *tmp;

    if (setting->bank)
        gabor_destroy(setting);

    i = _gabor_get_optimalsz(img_w);
    j = _gabor_get_optimalsz(img_h);
    opt_sz = i > j ? i : j;
    assert( opt_sz != 0 );
	
    setting->bank = fftw_malloc(sizeof(struct gabor_bank));
    assert( setting->bank != NULL );

    pbank = (struct gabor_bank*)setting->bank;
    pbank->opt_sz = opt_sz;
	real_num = setting->orientation_num * setting->scale_num;
    imag_num = real_num;
    len = opt_sz * opt_sz;
    sz  = sizeof(fftw_complex*) * (real_num + imag_num);
    sz += sizeof(fftw_complex) * (real_num + imag_num + 3) * len;
    sz += sizeof(double) * len * 2;

    pbank->pmem = fftw_malloc(sz);
    assert( pbank->pmem != NULL);

    memset(pbank->pmem, 0, sz);

    pbank->filters = (fftw_complex**)pbank->pmem;
    tmp = (fftw_complex*)(pbank->filters + real_num + imag_num);
    for (i = 0; i < real_num + imag_num; i += 2){
        /* put real part at [0], and imag part at [1] */
        *(pbank->filters + i + 0) = tmp; tmp += len;
        *(pbank->filters + i + 1) = tmp; tmp += len;
    }

    pbank->img_dfted = tmp; tmp += len;
    pbank->real_mult = tmp; tmp += len;
    pbank->imag_mult = tmp; tmp += len;

    pbank->buff0 = (double*)tmp;
    pbank->buff1 = (double*)(pbank->buff0 + len);

    tmp = NULL;

    pbank->plan_forward = fftw_plan_dft_r2c_2d(opt_sz,
                                               opt_sz,
                                               pbank->buff0,
                                               pbank->img_dfted,
                                               FFTW_MEASURE);

    pbank->plan_backward = fftw_plan_dft_c2r_2d(opt_sz,
                                                opt_sz,
                                                pbank->real_mult,
                                                pbank->buff0,
                                                FFTW_MEASURE);

    real_kernel_buff = (double*)fftw_malloc( sizeof(double) * len );
	imag_kernel_buff = (double*)fftw_malloc( sizeof(double) * len );
    assert( real_kernel_buff != NULL && 
		imag_kernel_buff != NULL );

    memset(real_kernel_buff, 0, sizeof(double) * len);
	memset(imag_kernel_buff, 0, sizeof(double) * len);

	k = 0;
    for (i = 0; i < setting->scale_num; ++i) {
        for (j = 0; j < setting->orientation_num; ++j) {
            _gabor_mk_kernel(setting->scales[ i ],
                             setting->orientations[ j ],
                             real_kernel_buff,
                             imag_kernel_buff,
                             opt_sz);

            _gabor_pad_kernel(real_kernel_buff, opt_sz, opt_sz,
                              pbank->buff0, opt_sz, opt_sz);

            _gabor_pad_kernel(imag_kernel_buff, opt_sz, opt_sz,
                              pbank->buff1, opt_sz, opt_sz);

            fftw_execute_dft_r2c(pbank->plan_forward,
                                 pbank->buff0,
                                 pbank->filters[ k++ ]);
            fftw_execute_dft_r2c(pbank->plan_forward,
                                 pbank->buff1,
                                 pbank->filters[ k++ ]);
        }
    }

    fftw_free( real_kernel_buff );
	fftw_free( imag_kernel_buff );
    return 0;
}

int gabor_extract(const struct gabor_setting *setting,
                  const unsigned char *image_data,
                  int w,
                  int ws,
                  int h,
                  double *feat_vec,
                  int feat_len) {

    struct gabor_bank *pbank = (struct gabor_bank*)setting->bank;
    double *buff0 = pbank->buff0;
    double *buff1 = pbank->buff1;
    double *pfeat = feat_vec;
    int opt_sz = pbank->opt_sz, len = pbank->opt_sz * pbank->opt_sz;
    int scale_num = setting->scale_num;
    int orientation_num = setting->orientation_num;
    fftw_complex *real_kernel, *imag_kernel;
    int i,j;

    /* set all buffer to 0 */
    memset(pbank->img_dfted,
           0,
           sizeof(fftw_complex) * (len * 3) + sizeof(double) * (len * 2));

    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            /*@TODO fast conversion */
            buff0[i*opt_sz + j] = image_data[i*ws + j];
        }
    }

    fftw_execute_dft_r2c(pbank->plan_forward,
                         buff0,
                         pbank->img_dfted);

    memset(buff0, 0, sizeof(double) * len);

    for (i = 0; i < scale_num; ++i) {
        for (j = 0; j < orientation_num; ++j) {

            /* @TODO if need to zero buff0 and buff1 every time ? */

            real_kernel = pbank->filters[ 2 * (i * orientation_num + j) + 0];
            imag_kernel = pbank->filters[ 2 * (i * orientation_num + j) + 1];
            /* multiply */
            _gabor_mult_scale(pbank->img_dfted,
                        real_kernel,
                        pbank->real_mult,
                        opt_sz);
            _gabor_mult_scale(pbank->img_dfted,
                        imag_kernel,
                        pbank->imag_mult,
                        opt_sz);

            /* inverse */
            fftw_execute_dft_c2r(pbank->plan_backward,
                                 pbank->real_mult, buff0);

            fftw_execute_dft_c2r(pbank->plan_backward,
                                 pbank->imag_mult, buff1);
            /* 'MAG' */
            _gabor_calc_mag(buff0, buff1, buff0, w, opt_sz, h);

            /* extract: custom extract method should guarantee,
               the feat_vec is long enough.
               !note: this function changes pfeat */
            _gabor_custom_extract(buff0, &pfeat, w, opt_sz, h);
        }
    }

    return 0;
}

void gabor_test(struct gabor_setting *setting, int fidx,
	const unsigned char *image_data, int w, int ws, int h, double *vec) {

		int i, j;
		fftw_complex *real_kernel, *imag_kernel;
		struct gabor_bank *pbank = (struct gabor_bank*)setting->bank;
		double *buff0 = pbank->buff0;
		double *buff1 = pbank->buff1;
		int opt_sz = pbank->opt_sz;
		int len = opt_sz * opt_sz;

		real_kernel = pbank->filters[2 * fidx + 0];
		imag_kernel = pbank->filters[2 * fidx + 1];

		for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            /*@TODO fast conversion */
            buff0[i*opt_sz + j] = image_data[i*ws + j];
        }
    }

    fftw_execute_dft_r2c(pbank->plan_forward,
                         buff0,
                         pbank->img_dfted);

    memset(buff0, 0, sizeof(double) * len);
	
	 _gabor_mult_scale(pbank->img_dfted,
                        real_kernel,
                        pbank->real_mult,
                        opt_sz);
	 _gabor_mult_scale(pbank->img_dfted,
		 imag_kernel,
		 pbank->imag_mult,
		 opt_sz);

	 /* inverse */
	 fftw_execute_dft_c2r(pbank->plan_backward,
		 pbank->real_mult, buff0);

	 fftw_execute_dft_c2r(pbank->plan_backward,
		 pbank->imag_mult, buff1);
	 /* 'MAG' */
	 _gabor_calc_mag(buff0, buff1, buff0, w, opt_sz, h);

	 for (i = 0; i < h; ++i) {
		 for (j = 0; j < w; ++j) {
			 vec[i*w + j] = buff0[i*opt_sz + j];
		 }
	 }
	 _gabor_minmax(vec, w*h);
}

EXTERN_END