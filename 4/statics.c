/**
 * Temporary put static utils here.
 */
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int _gabor_get_optimalsz(int n) {
    int m, newsz = n;
    for (;;) {
        m = newsz;
        while ((m % 2) == 0) m /= 2;
        while ((m % 3) == 0) m /= 3;
        while ((m % 5) == 0) m /= 5;
        while ((m % 7) == 0) m /= 7;
        if (m <= 1) break;
        newsz ++;
    }
    return newsz;
}

static void _gabor_scale(double *vec, int n, double scale) {
    double *p = vec, *end = p + n;
    do {
        *p ++ *= scale;
    } while (p != end);
}


#ifndef GABOR_PI
#define GABOR_PI 3.1415926535897932384626433832795
#endif
static void _gabor_mk_kernel(int iNu, int iMu,
                             double *real,
                             double *imag,
                             int opt_sz) {
	int x, y,i,j;
    double dReal;
    double dImag;
    double dTemp1, dTemp2, dTemp3;
	double Sigma = 2 * GABOR_PI;
	double Kmax = GABOR_PI / 2;
	double dFreq = sqrt(2.0);
	double K = Kmax / pow(dFreq, (double)iNu);
	double Phi = iMu * (GABOR_PI/8);
    for (i = 0; i < opt_sz; i++)
    {
		y = i-(opt_sz-1)/2;
        for (j = 0; j < opt_sz; j++)
        {
            x = j-(opt_sz-1)/2;
            dTemp1 = (pow(K,2)/pow(Sigma,2))*exp(-(pow((double)x,2)+pow((double)y,2))*pow(K,2)/(2*pow(Sigma,2)));
            dTemp2 = cos(K*cos(Phi)*y + K*sin(Phi)*x) - exp(-(pow(Sigma,2)/2));
            dTemp3 = sin(K*cos(Phi)*y + K*sin(Phi)*x);
            dReal = dTemp1*dTemp2;
            dImag = dTemp1*dTemp3; 
 
			real[i*opt_sz + j] = dReal;
			imag[i*opt_sz + j] = dImag;
        } 
    }
}


static void _gabor_pad_kernel(double *src, int sw, int sh,
                              double *dst, int dw, int dh) {
    int si,sj,di,dj;

    int cent_x = (sw - 1) / 2, cent_y = (sh - 1) / 2;

    double *psrc = src;
    double *pdst = dst;

    /* kernel size should not bigger than buffer size */
    if (sw > dw || sh > dh) return ;

    for (di = 0, si = cent_y; si < sh; ++si, ++di)
        for (dj = 0, sj = cent_x; sj < sw; ++sj, ++dj)
            pdst[di*dw + dj] = psrc[si*sw + sj];

    for (di = dh-1, si = cent_y - 1; si >= 0; --si, --di)
        for (dj = 0, sj = cent_x + 1; sj < sw; ++sj, ++dj)
            pdst[di*dw + dj] = psrc[si*sw + sj];

    for (di = 0, si = cent_y + 1; si < sh; ++si, ++di)
        for (dj = dw - 1, sj = cent_x - 1; sj >= 0; --sj, --dj)
            pdst[di*dw + dj] = psrc[si*sw + sj];

    for (di = dh - 1, si = cent_y - 1; si >= 0; --si, --di)
        for (dj = dw - 1, sj = cent_x - 1; sj >= 0; --sj, --dj)
            pdst[di*dw + dj] = psrc[si*sw + sj];
}

static void _gabor_mult_scale(const fftw_complex *src0,
                        const fftw_complex *src1,
                        fftw_complex *dst,
                        int opt_sz) {
    int i, j, ij, half = opt_sz/2 + 1;
    const fftw_complex *psrc0 = src0, *psrc1 = src1;
    fftw_complex *pdst = dst;
    double k1, k2, k3,a,b,c,d;
	double scale = 1.0 / (opt_sz * opt_sz);

    for (i = 0; i < opt_sz; ++i) {
		for (j = 0; j < half; ++j) {
			ij = i * half + j;

			a = psrc0[ij][0];
			b = psrc0[ij][1];
			c = psrc1[ij][0];
			d = psrc1[ij][1];

			k1 = a * (c + d);
			k2 = d * (a + b);
			k3 = c * (b - a);

			pdst[ij][0] = (k1 - k2) * scale;
			pdst[ij][1] = (k1 + k3) * scale ;
		}
	}
}

static void _gabor_calc_mag(const double *real,
                            const double *imag,
                            double *dst,
                            int w, int ws, int h) {
    int i, j, idx;
    const double
        *preal = real,
        *pimag = imag;
    double *pdst = dst;

    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            idx = i*ws +j;
            pdst[ idx ] = sqrt(preal[idx] * preal[idx] +
                               pimag[idx] * pimag[idx]);
        }
    }
}

/**
 * !note: this minmax should not perform on
 * a mid-result like buff0 or buff1, 'cause
 * the valid part is probably in the top-left
 * part.
 */
static void _gabor_minmax(double *vec, int n) {
    double min,max, t, *pvec = vec, *end = vec + n;
    min = max = *pvec;
    do {
        t = *pvec++;
        min = t > min ? min : t;
        max = t < max ? max : t;
    } while (pvec != end);

    pvec = vec;
    max = max - min;
    if (max != 0) max = 1.0 / max;
    do {
        *pvec -= min;
        *pvec *= max;
        pvec ++;
    } while (pvec != end);
}

static void _gabor_minmax_rect(double *vec, int w, int ws, int h) {
    int i, j, idx;
    double t, min, max;
    min = max = *vec;
    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            idx = i * ws + j;
            t = vec[ idx ];
            min = t > min ? min : t;
            max = t < max ? max : t;
        }
    }

    max = max - min;
    if (max != 0) {
        max = 1.0/max;
        for (i = 0; i < h; ++i) {
            for (j = 0; j < w; ++j) {
                idx = i * ws + j;
                vec[ idx ] -= min;
                vec[ idx ] *= max;
            }
        }
    }
}

/**
 * You probably should use your own routine to
 * extract dense features, here is a simple
 * example.
 */
static void _gabor_custom_extract(double *src,
                                  double **dst,
                                  int w,
                                  int ws,
                                  int h) {
    int i,j;
    int step_x = 10, step_y = 10;
    double *pdst = *dst;

    _gabor_minmax_rect( src , w, ws, h);

    for (i = step_y; i < h; i += step_y) {
        for (j = step_x; j < w; j += step_x) {
            *pdst++ = src[i*ws + j];
        }
    }
}
