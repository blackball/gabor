/**
 * Make gabor kernel.
 *
 * @blackball
 */

/**
 * make kernel
 * @kw kernel width
 * @kh kernel height
 * @
 */
#define GABOR_PI 3.1415926535897932384626433832795
void gabor_mkKernel(int kw,int kh,
                    double sigma, /* default = 1*/
                    double gamma,
                    double theta,
                    double lambda,
                    double psi,
                    double *real,
                    double *imag) {
  /* half kernel width and height, should be odd, if it's not odd,
   then get the biggest odd sub rect, and fill the border with 0s */
  int hkw = kw / 2, hkh = kh / 2;
  
  double _theta = theta * (GABOR_PI/180);
  double _psi = psi * (GABOR_PI/180);
  double _gamma = gamma * gamma;
  
  /* what 2.0 means ?! */
  double x_, y_;
  double e_;

  /* normalize */
  double xdelta_ = 2.0 / hkw; 
  double ydelta_ = 2.0 / hkh;
  
  for (int iy = -hkh; iy < hkh; ++iy)
    for (int ix = -hkw; ix < hkw; ++ix) {
      /**
       * @TODO much optimization could be performed!
       */
      x_ = ix * xdelta_ * cos(_theta)  + iy * ydelta_ * sin(_theta);
      y_ = -ix * xdelta_ * sin(_theta) + iy * ydelta_ * con(_theta);
      e_ = exp(- ((x_* x_) + _gamma * (y_ * y_)) / (2*GABOR_PI));

      
      real(hkw + ix, hkh + iy) = e_ * cos(2*(GABOR_PI) * x_/lambda + _phi);
      imag(hkw + ix, hkh + iy) = e_ * sin(2*(GABOR_PI) * x_/lambda + _phi);
    } 
}

void gabor_init() {
  
}
