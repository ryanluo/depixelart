#ifndef IP_H
#define IP_H


#include "common.h"
#include "image.h"


/*
* IMPORTANT - DO NOT CHANGE THE INTERFACES DEFINED HERE - IMPORTANT
*
* It's very important that the interface to your program does not
* change from what is provided, so that automated testing scripts
* will work.  If you add additional requests for input, these scripts
* will no longer work and it will negatively affect your grade.  Each
* method has sufficient information to carry out the required task.
*
* The only method you are may get more user input in is ip_warp().
* This option will not be tested automatically, since everyone's
* will be different, so you may ask for whatever input is necessary.
* To support multiple warps, have ip_warp() print a menu of options.
*
* Of course, you may add whatever other functions you may need.
*/


Image*	ip_blur_box (Image* src, int size);
Image*	ip_blur_gaussian (Image* src, int size, double sigma);
Image*	ip_blur_triangle (Image* src, int size);
Image*	ip_brighten (Image* src, double alpha);
Image*	ip_color_shift (Image* src);
Image*	ip_composite (Image* src1, Image* src2, Image* mask);
Image*	ip_contrast (Image* src, double alpha);
Image*  ip_convolve(Image* src, int size, double* kernel);
Image*	ip_crop (Image* src, int x0, int y0, int x1, int y1);
Image*	ip_edge_detect (Image* src);
Image*	ip_extract (Image* src, int channel);
Image*	ip_fun_warp (Image* src, int samplingMode);
Image*	ip_grey (Image* src);
Image*  ip_image_shift(Image* src, double dx, double dy);
Image*  ip_interpolate(Image* src1, Image* src2, double alpha);
Image*	ip_invert (Image* src);
Image*	ip_misc(Image* src, ImageGraph* imageGraph, vector<vector<bool>> * similarityGraph, const int blockSizeX, const int blockSizeY, const int pixelSize);
Image*	ip_quantize_simple (Image* src, int bitsPerChannel);
Image*	ip_quantize_ordered (Image* src, int bitsPerChannel);
Image*	ip_quantize_fs (Image* src, int bitsPerChannel);
Pixel	ip_resample_nearest(Image* src, double x, double y);
Pixel	ip_resample_bilinear(Image* src, double x, double y);
Pixel   ip_resample_gaussian(Image* src, double x, double y, int filtersize, double sigma);
Image*	ip_rotate (Image* src, double theta, int x, int y, int samplingMode, int gaussianFilterSize, double gaussianSigma);
Image*	ip_saturate (Image* src, double alpha);
Image*	ip_scale (Image* src, double x, double y, int samplingMode, int gaussianFilterSize, double gauusianSigma);
Image*	ip_threshold (Image* src, double cutoff);




#endif // IP_H
