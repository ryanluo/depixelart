#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
 *By
 *Zhang Zehao(Michael)
 *Zhu Huangjian(Sean)
 */

/*
* convolve with a box filter
*/
Image* ip_blur_box (Image* src, int size)
{
    int area = sqr(size);
    double kernel[area];
    for(int i=0;i<area;i++){
        kernel[i] = 1.0/area;
    }
    return ip_convolve(src, size, kernel);
}


/*
* convolve with a gaussian filter
*/
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    int start = (size-1.0)/2.0;
    int area = sqr(size);
    double kernel[area];
    double sum = 0.0;
    int x, y;
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            x = j-start;
            y = i-start;
            kernel[i*size+j] = 1.0/pow(M_E,(sqr(x)+sqr(y))/(2*sqr(sigma)));
            sum = sum+kernel[i*size+j];
        }
    }
    for(int i=0; i<area; i++){
        kernel[i] = kernel[i]/sum;
    }
    return ip_convolve(src, size, kernel);
}


/*
* convolve with a triangle filter
*/
Image* ip_blur_triangle (Image* src, int size)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* interpolate with a black image
*/
Image* ip_brighten (Image* src, double alpha)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    Image* black = new Image(wd,ht,bt);
    return ip_interpolate(black, src, alpha);
}


/*
* shift colors
*/
Image* ip_color_shift(Image* src)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    Image* shift = new Image(wd,ht,bt);
    double srcR,srcG,srcB;
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            srcR = src->getPixel(i, j, 0);
            srcG = src->getPixel(i, j, 1);
            srcB = src->getPixel(i, j, 2);
            shift->setPixel(i, j, 1, srcR);
            shift->setPixel(i, j, 2, srcG);
            shift->setPixel(i, j, 0, srcB);
        }
    }
	return shift;
}


/*
* use a mask image for a per-pixel alpha value to perform
* interpolation with a second image
*/
Image* ip_composite (Image* src1, Image* src2, 
					 Image* mask)
{
    int wd = src1->getWidth();
    int ht = src1->getHeight();
    int bt = src1->getBits();
    double src1R, src1G,src1B,src2R,src2G,src2B,maskValue,compR,compG,compB;
    Image* comp = new Image(wd,ht,bt);
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            maskValue = mask->getPixel(i,j,0);
            src1R = src1->getPixel(i, j, 0);
            src1G = src1->getPixel(i, j, 1);
            src1B = src1->getPixel(i, j, 2);
            src2R = src2->getPixel(i, j, 0);
            src2G = src2->getPixel(i, j, 1);
            src2B = src2->getPixel(i, j, 2);
            compR = (1-maskValue)*src1R + maskValue*src2R;
            compG = (1-maskValue)*src1G + maskValue*src2G;
            compB = (1-maskValue)*src1B + maskValue*src2B;
            compR = clamp(compR, 0, 1);
            compG = clamp(compG, 0, 1);
            compB = clamp(compB, 0, 1);
            comp->setPixel(i, j, 0, compR);
            comp->setPixel(i, j, 1, compG);
            comp->setPixel(i, j, 2, compB);
        }
    }
    return comp;
}


/*
* interpolate with the average intensity of the src image
*/
Image* ip_contrast (Image* src, double alpha)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    int amountPixels = wd*ht;
    double cR = 0.2126;
    double cG = 0.7152;
    double cB = 0.0722;
    double srcRTotal=0, srcGTotal=0, srcBTotal=0, srcRAv,srcGAv,srcBAv;
    Image* grey = new Image(wd,ht,bt);
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            srcRTotal = srcRTotal+src->getPixel(i, j, 0);
            srcGTotal = srcGTotal+src->getPixel(i, j, 1);
            srcBTotal = srcBTotal+src->getPixel(i, j, 2);
        }
    }
    srcRAv = srcRTotal/amountPixels;
    srcGAv = srcGTotal/amountPixels;
    srcBAv = srcGTotal/amountPixels;
    double gValue = cR*srcRAv+cG*srcGAv+cB*srcBAv;
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            grey->setPixel(i, j, 0, gValue);
            grey->setPixel(i, j, 1, gValue);
            grey->setPixel(i, j, 2, gValue);
        }
    }
    return ip_interpolate(src, grey, alpha);

}



/*
* convolve an image with a kernel
*/
Image* ip_convolve (Image* src, int size, double* kernel )
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    Image* conv = new Image(wd,ht,bt);
    double sumR=0.0,sumG=0.0,sumB=0.0;
    double srcR,srcG,srcB;
    double weight;
    int startX,startY,pX,pY;
    for(int i=0;i<wd;i++){
        startX = i-((size-1)/2);
        for(int j=0;j<ht;j++){
            startY = j-((size-1)/2);
            for(int a=0;a<size;a++){
                for(int b=0;b<size;b++){
                    pX = startX+a;
                    pY = startY+b;
                    if(pX>=0&&pY>=0&&pX<wd&&pY<ht){
                        weight = *(kernel+(b*size+a));
                        srcR = src->getPixel(pX, pY, 0);
                        srcG = src->getPixel(pX, pY, 1);
                        srcB = src->getPixel(pX, pY, 2);
                        sumR = sumR + srcR*weight;
                        sumG = sumG + srcG*weight;
                        sumB = sumB + srcB*weight;
                    }
                }
            }
            sumR = clamp(sumR,0.0,1.0);
            sumG = clamp(sumG,0.0,1.0);
            sumB = clamp(sumB,0.0,1.0);
            conv->setPixel(i, j, 0, sumR);
            conv->setPixel(i, j, 1, sumG);
            conv->setPixel(i, j, 2, sumB);
            sumR=0.0;
            sumG=0.0;
            sumB=0.0;
        }
    }
    return conv;
}



/*
*  create cropped version of image
*/
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    if (x1 < src->getWidth() && y1 < src ->getHeight() && x0 >=0 && y0 >= 0){
        int wd = x1 - x0;
        int ht = y1 - y0;
        int bt = src -> getBits();
        Pixel sample;
        Image *newImage = new Image(wd,ht,bt);
        for (int i = 0; i < wd; i++){
            for (int j = 0; j < ht; j++){
                sample = ip_resample_nearest(src, x0+i, y0+j);
                newImage -> setPixel(i, j, sample);
            }
        }
        return newImage;
    }else{
        return NULL;
    }
    
}
/*
* convolve with an edge detection kernel
*/
Image* ip_edge_detect (Image* src)
{
    int size = 3;
    int area = sqr(size);
    int mid = 4;
    double kernel[9];
    for(int i=0;i<area;i++){
        if(i==mid){
            kernel[i] = 8;
        }else{
            kernel[i] = -1;
        }
    }
    return ip_convolve(src, size, kernel);
}


/*
* extract channel of input image
*/
Image* ip_extract (Image* src, int channel)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    double srcC;
    Image* extract = new Image(wd,ht,bt);
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            srcC = src->getPixel(i, j, channel);
            extract->setPixel(i, j, channel, srcC);
        }
    }
    return extract;
}


/*
* create your own fun warp
*/
Image* ip_fun_warp (Image* src, int samplingMode)
{
    double magnitude;
    double frequency;
    cout<<"Please input the magnitude of the wave[0,1) "<<endl;
    cin>>magnitude;
    cout<<"Please input the frequency of the wave "<<endl;
    cin>>frequency;
    magnitude = 1 / (1 - magnitude);
    frequency = frequency;
    int wd = src -> getWidth();
    int ht = src -> getHeight();
    int bt = src -> getBits();
    double newX;
    Pixel pix;
    Image *fun = new Image(wd,ht,bt);
    for(int i = 0; i < wd; i++){
        for(int j=0;j<ht;j++){
            newX = i - (magnitude * sin(frequency * (double)j));
            if(newX > wd-1){
                newX = newX - (wd-1);
            }else if(newX < 0){
                newX = newX + (wd - 1);
            }
            if(samplingMode == 0){
                pix = ip_resample_nearest(src, newX, j);
            }else if(samplingMode == 1){
                pix = ip_resample_bilinear(src, newX, j);
            }else{
                pix = ip_resample_gaussian(src, newX, j, 4, 1);
            }
            fun -> setPixel(i, j, pix);
        }
    }
	//  ask user for input parameters here including resampling method and, 
	//  if gaussian resampling is used, its filtersize and sigma
	//  if you implement more than one warp, you should ask the 
	//  user to chose the one to perform here too!
    
    

	return fun;
}
/*
* create a new image with values equal to the psychosomatic intensities
* of the source image
*/
Image* ip_grey (Image* src)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    double cR = 0.2126;
    double cG = 0.7152;
    double cB = 0.0722;
    double srcR,srcG,srcB,gValue;
    Image* grey = new Image(wd,ht,bt);
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            srcR = src->getPixel(i, j, 0);
            srcG = src->getPixel(i, j, 1);
            srcB = src->getPixel(i, j, 2);
            gValue = cR*srcR+cG*srcG+cB*srcB;
            grey->setPixel(i, j, 0, gValue);
            grey->setPixel(i, j, 1, gValue);
            grey->setPixel(i, j, 2, gValue);
        }
    }
    return grey;
}


/*
*  shift image by dx and dy (modulo width & height)
*/

Image* ip_image_shift (Image* src, double dx, double dy)
{
    int wd = src -> getWidth();
    int ht = src -> getHeight();
    int bt = src -> getBits();
    Image *shift = new Image(wd,ht,bt);
    Pixel pix;
    int modX, modY;
    for(int i = 0; i < wd; i++){
        for(int j = 0; j < ht; j++){
            pix = src -> getPixel(i, j);
            modX = (i + (int)(dx + 0.5)) % wd;
            modY = (j + (int)(dy + 0.5)) % ht;
            shift -> setPixel(modX, modY, pix);
        }
    }
	return shift;
}
/*
* interpolate an image with another image
*/
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    //because the two images have the same dimensions and bits per pixel by defualt
    //we just take the width,length, and BPP from the first image.
    int wd = src1->getWidth();
    int ht = src1->getHeight();
    int bt = src1->getBits();
    Image* inter = new Image(wd,ht,bt);
    double src1R,src1G,src1B,src2R,src2G,src2B,interR,interG,interB;
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            src1R = src1->getPixel(i, j, 0);
            src1G = src1->getPixel(i, j, 1);
            src1B = src1->getPixel(i, j, 2);
            src2R = src2->getPixel(i, j, 0);
            src2G = src2->getPixel(i, j, 1);
            src2B = src2->getPixel(i, j, 2);
            interR = alpha*src1R + (1-alpha)*src2R;
            interG = alpha*src1G + (1-alpha)*src2G;
            interB = alpha*src1B + (1-alpha)*src2B;
            interR = clamp(interR, 0, 1);
            interG = clamp(interG, 0, 1);
            interB = clamp(interB, 0, 1);
            inter->setPixel(i, j, 0, interR);
            inter->setPixel(i, j, 1, interG);
            inter->setPixel(i, j, 2, interB);
        }
    }
    return inter;
}
/*
* invert input image
*/
Image* ip_invert (Image* src)
{
	/*to make the program run more efficient, we
     *replaced the values of the intensity of the
     *second image and of the alpha directly into
     *the computation of interpolation*/
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    Image* invert = new Image(wd,ht,bt);
    double srcR,srcG,srcB,invR,invG,invB;
    double grey = 0.5,alpha=-1;
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            srcR = src->getPixel(i, j, 0);
            srcG = src->getPixel(i, j, 1);
            srcB = src->getPixel(i, j, 2);
            invR = alpha*srcR + (1-alpha)*grey;
            invG = alpha*srcG + (1-alpha)*grey;
            invB = alpha*srcB + (1-alpha)*grey;
            invR = clamp(invR, 0, 1);
            invG = clamp(invG, 0, 1);
            invB = clamp(invB, 0, 1);
            invert->setPixel(i, j, 0, invR);
            invert->setPixel(i, j, 1, invG);
            invert->setPixel(i, j, 2, invB);
        }
    }
    return invert;

}


/*
* define your own filter
* you need to request any input paraters here, not in control.cpp
*/

Image* ip_misc(Image* src)
{
    int choice;
    cout<<"Please select effect: "<<endl;
    cout<<"      1.gamma correction"<<endl;
    cout<<"      2.bilateral filter"<<endl;
    cout<<"      3.median filter"<<endl;
    cout<<"      4.sobel operator"<<endl;
    cout<<"      5.video"<<endl;
    cin>>choice;
    //gamma correction
    if(choice == 1){
        double gammaV;
        cout<<"Please input gamma value(>0): "<<endl;
        cin >> gammaV;
        int wd = src -> getWidth();
        int ht = src -> getHeight();
        int bt = src -> getBits();
        Image *gamma = new Image(wd,ht,bt);
        double srcR, srcG, srcB, corR, corG, corB;
        for(int i=0;i<wd;i++){
            for(int j=0;j<ht;j++){
                srcR = src -> getPixel(i, j, 0);
                srcG = src -> getPixel(i, j, 1);
                srcB = src -> getPixel(i, j, 2);
                corR = pow(srcR, gammaV);
                corG = pow(srcG, gammaV);
                corB = pow(srcB, gammaV);
                gamma -> setPixel(i, j, 0, corR);
                gamma -> setPixel(i, j, 1, corG);
                gamma -> setPixel(i, j, 2, corB);
            }
        }
        return gamma;
    }
    //bilateral filter
    else if(choice == 2){
        int size;
        double sigma1;
        double sigma2;
        cout<<"Please input filter size(odd): ";
        cin>>size;
        cout<<"Please input sigma for intensity: ";
        cin>>sigma1;
        cout<<"Please input sigma for distance: ";
        cin>>sigma2;
        int wd = src -> getWidth();
        int ht = src -> getHeight();
        int bt = src -> getBits();
        Image *bilateral = new Image(wd,ht,bt);
        int startX, startY;
        double srcR, srcG, srcB, centerR, centerG, centerB,  weights, weightiR, weightiG, weightiB, finalR, finalG,finalB;
        for (int x = 0; x < wd; x++){
            for (int y = 0; y < ht; y++){
                double sumR=0, sumG=0, sumB=0, sumWeightR=0, sumWeightG=0, sumWeightB=0;
                startX = x - (size/2);
                startY = y - (size/2);
                centerR = src -> getPixel(x, y, 0);
                centerG = src -> getPixel(x, y, 1);
                centerB = src -> getPixel(x, y, 2);
                for (int i = 0; i < size; i++){
                    for (int j = 0; j < size; j++){
                        if(startX + i >= 0 && startX + i < wd && startY + j >=0 && startY + j < ht){
                            srcR = src -> getPixel(startX + i, startY + j, 0);
                            srcG = src -> getPixel(startX + i, startY + j, 1);
                            srcB = src -> getPixel(startX + i, startY + j, 2);
                            weights = (1.0/pow(M_E,(sqr(x - (startX + i))+sqr(y - (startY + j)))/(2*sqr(sigma2)))) ;
                            //cout<<"sigma1"<<sigma1<<endl;
                            weightiR = (1.0/pow(M_E,(sqr(centerR - srcR))/(2*sqr(sigma1))));
                            sumWeightR += weightiR*weights;
                            sumR += weightiR*weights*srcR;
                            weightiG = (1.0/pow(M_E,(sqr(centerG - srcG))/(2*sqr(sigma1))));
                            sumWeightG += weightiG*weights;
                            sumG += weightiG*weights*srcG;
                            weightiB = (1.0/pow(M_E,(sqr(centerB - srcB))/(2*sqr(sigma1))));
                            sumWeightB += weightiB*weights;
                            sumB += weightiB*weights*srcB;
                        }
                    }
                }
                finalR = sumR / sumWeightR;
                finalG = sumG / sumWeightG;
                finalB = sumB / sumWeightB;
                finalR = clamp(finalR, 0, 1);
                finalG = clamp(finalG, 0, 1);
                finalB = clamp(finalB, 0, 1);
                bilateral -> setPixel(x, y, 0, finalR);
                bilateral -> setPixel(x, y, 1, finalG);
                bilateral -> setPixel(x, y, 2, finalB);
            }
        }
        return bilateral;
        
    }
    //median filter
    else if(choice == 3){
        int size;
        cout<<"Please input filter size(odd): ";
        cin>>size;
        int wd = src -> getWidth();
        int ht = src -> getHeight();
        int bt = src -> getBits();
        Image *median = new Image(wd,ht,bt);
        int startX, startY;
        double rValues[size*size];
        double gValues[size*size];
        double bValues[size*size];
        double srcR, srcG, srcB;
        int index;
        for(int x=0; x<wd;x++){
            for(int y=0; y<ht;y++){
                startX = x -(size/2);
                startY = y -(size/2);
                index = 0;
                for(int i=0; i<size; i++){
                    for(int j=0; j<size; j++){
                        if(startX + i >= 0 && startX + i < wd && startY + j >=0 && startY + j < ht){
                            srcR = src -> getPixel(startX + i, startY + j, 0);
                            srcG = src -> getPixel(startX + i, startY + j, 1);
                            srcB = src -> getPixel(startX + i, startY + j, 2);
                            rValues[index] = srcR;
                            gValues[index] = srcG;
                            bValues[index] = srcB;
                            index++;
                        }
                    }
                }
                sort(rValues, rValues+index+1);
                sort(gValues, gValues+index+1);
                sort(bValues, bValues+index+1);
                median -> setPixel(x, y, 0, rValues[(index+1)/2]);
                median -> setPixel(x, y, 1, gValues[(index+1)/2]);
                median -> setPixel(x, y, 2, bValues[(index+1)/2]);
            }
        }
        return median;
    }
    //sobel operator
    else if(choice == 4){
        int wd = src -> getWidth();
        int ht = src -> getHeight();
        int bt = src -> getBits();
        int size =3;
        double kernelv[] = {-1,0,1,-2,0,2,-1,0,1};
        double kernelh[] = {-1,-2,-1,0,0,0,1,2,1};
        Image *sobel = new Image(wd,ht,bt);
        double sumRh=0.0,sumGh=0.0,sumBh=0.0, sumRv=0.0,sumGv=0.0,sumBv=0.0,finalR,finalG,finalB;
        double srcR,srcG,srcB;
        double weighth,weightv;
        int startX,startY,pX,pY;
        for(int i=0;i<wd;i++){
            startX = i-((size-1)/2);
            for(int j=0;j<ht;j++){
                startY = j-((size-1)/2);
                for(int a=0;a<size;a++){
                    for(int b=0;b<size;b++){
                        pX = startX+a;
                        pY = startY+b;
                        if(pX>=0&&pY>=0&&pX<wd&&pY<ht){
                            weighth = *(kernelh+(b*size+a));
                            weightv = *(kernelv+(b*size+a));
                            srcR = src->getPixel(pX, pY, 0);
                            srcG = src->getPixel(pX, pY, 1);
                            srcB = src->getPixel(pX, pY, 2);
                            sumRh = sumRh + srcR*weighth;
                            sumGh = sumGh + srcG*weighth;
                            sumBh = sumBh + srcB*weighth;
                            sumRv = sumRv + srcR*weightv;
                            sumGv = sumGv + srcG*weightv;
                            sumBv = sumBv + srcB*weightv;
                        }
                    }
                }
                finalR = sqrt(sqr(sumRh)+sqr(sumRv));
                finalG = sqrt(sqr(sumGh)+sqr(sumGv));
                finalB = sqrt(sqr(sumBh)+sqr(sumBv));
                finalR = clamp(finalR,0,1);
                finalG = clamp(finalG,0,1);
                finalB = clamp(finalB,0,1);
                sobel -> setPixel(i, j, 0, finalR);
                sobel -> setPixel(i, j, 1, finalG);
                sobel -> setPixel(i, j, 2, finalB);
                sumRh=0.0;
                sumGh=0.0;
                sumBh=0.0;
                sumRv=0.0;
                sumGv=0.0;
                sumBv=0.0;
            }
        }
        return sobel;
    }
    //video filter
    else if(choice == 5){
        int bandWidth;
        cout<<"Please input bandwidth(integers > 0): "<<endl;
        cin>>bandWidth;
        int wd = src -> getWidth();
        int ht = src -> getHeight();
        int bt = src -> getBits();
        double srcR, srcG, srcB;
        src = ip_brighten(src, -0.3);
        src = ip_saturate(src, 1.2);
        Image *video = new Image(wd,ht,bt);
        for(int i=0; i<wd; ++i){
            for (int j=0; j<ht; j+=bandWidth){
                if(j%(3*bandWidth) == 0){
                    for(int a=0; a<bandWidth; ++a){
                        if(j+a < ht){
                            srcR = src -> getPixel(i, j+a, 0);
                            srcG = src -> getPixel(i, j+a, 1);
                            srcB = src -> getPixel(i, j+a, 2);
                            video -> setPixel(i, j+a, 0, srcR);
                            video -> setPixel(i, j+a, 1, srcG*.6);
                            video -> setPixel(i, j+a, 2, srcB*.6);
                        }
                    }
                }else if(j%(3*bandWidth) == bandWidth){
                    for(int a=0; a<bandWidth; ++a){
                        if(j+a < ht){
                            srcR = src -> getPixel(i, j+a, 0);
                            srcG = src -> getPixel(i, j+a, 1);
                            srcB = src -> getPixel(i, j+a, 2);
                            video -> setPixel(i, j+a, 1, srcG);
                            video -> setPixel(i, j+a, 0, srcR*.6);
                            video -> setPixel(i, j+a, 2, srcB*.6);
                        }
                    }
                }else if(j%(3*bandWidth) == 2*bandWidth){
                    for(int a=0; a<bandWidth; ++a){
                        if(j+a < ht){
                            srcR = src -> getPixel(i, j+a, 0);
                            srcG = src -> getPixel(i, j+a, 1);
                            srcB = src -> getPixel(i, j+a, 2);
                            video -> setPixel(i, j+a, 2, srcB);
                            video -> setPixel(i, j+a, 0, srcR*.6);
                            video -> setPixel(i, j+a, 1, srcG*.6);
                        }
                    }
                }
            }
        }
        return video;
    }
	return NULL;
}




/*
* round each pixel to the nearest value in the new number of bits
*/
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    double srcR,srcG,srcB;
    Image* quan = new Image(wd,ht,bitsPerChannel);
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            srcR = src->getPixel(i, j, 0);
            srcG = src->getPixel(i, j, 1);
            srcB = src->getPixel(i, j, 2);
            quan->setPixel(i, j, 0, srcR);
            quan->setPixel(i, j, 1, srcG);
            quan->setPixel(i, j, 2, srcB);
        }
    }
    return quan;
}


/*
* dither each pixel to the nearest value in the new number of bits
* using a static 4x4 matrix
*/
Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    Image* quanOrd = new Image(wd,ht,bt);
    double srcR,srcG,srcB,offSet;
    int yOffset[4] = {1,0,1,0};
    int xOffset[4] = {0,1,1,0};
    int outPutBlockNum = (pow(2,bitsPerChannel) - 1) * 8;
    for(int x = 0; x < wd-1; x+=2){
        for(int y = 0; y < ht-1; y+=2){
            for(int i=0; i < 4; ++i){
                int xPos = x + xOffset[i];
                int yPos = y + yOffset[i];
                srcR = src->getPixel(xPos, yPos, 0);
                srcG = src->getPixel(xPos, yPos, 1);
                srcB = src->getPixel(xPos, yPos, 2);
                offSet = (3.0 - 2.0*i) / outPutBlockNum;
                srcR += offSet;
                srcG += offSet;
                srcB += offSet;
                srcR = clamp(srcR, 0, 1);
                srcG = clamp(srcG, 0, 1);
                srcB = clamp(srcB, 0, 1);
                quanOrd->setPixel(xPos, yPos, 0, srcR);
                quanOrd->setPixel(xPos, yPos, 1, srcG);
                quanOrd->setPixel(xPos, yPos, 2, srcB);
            }
        }
    }
    return ip_quantize_simple(quanOrd, bitsPerChannel);
}


/*
* dither each pixel to the nearest value in the new number of bits
* using error diffusion
*/
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    double Rerr[ht][wd];
    double Gerr[ht][wd];
    double Berr[ht][wd];
    for(int i=0;i<ht;i++){
        for(int j=0;j<wd;j++){
            Rerr[i][j] = 0.0;
            Gerr[i][j] = 0.0;
            Berr[i][j] = 0.0;
        }
    }
    double srcR,srcG,srcB,Rcopy,Gcopy,Bcopy, Rerror, Gerror, Berror;
    double levels = pow(2,bitsPerChannel);
    double b = 3.0/16.0, a = 7.0/16.0, x = 5.0/16.0, y = 1.0/16.0;
    Image* quan = new Image(wd,ht,bitsPerChannel);
    for(int j=0;j<ht;j++){
        for(int i=0;i<wd;i++){
            
            srcR = src->getPixel(i, j, 0);
            srcG = src->getPixel(i, j, 1);
            srcB = src->getPixel(i, j, 2);
            srcR = srcR + Rerr[j][i];
            srcG = srcG + Gerr[j][i];
            srcB = srcB + Berr[j][i];
            Rcopy = ((double)(int)(srcR * (levels - 1) + 0.5)) / (levels - 1.0);
            Gcopy = ((double)(int)(srcG * (levels - 1) + 0.5)) / (levels - 1.0);
            Bcopy = ((double)(int)(srcB * (levels - 1) + 0.5)) / (levels - 1.0);
            Rerror = srcR - Rcopy;
            Gerror = srcG - Gcopy;
            Berror = srcB - Bcopy;
            Rcopy = clamp(Rcopy,0,1);
            Gcopy = clamp(Gcopy,0,1);
            Bcopy = clamp(Bcopy,0,1);
            if(i > 0 && j < ht-1){
                Rerr[j+1][i-1] += b * Rerror;
                Gerr[j+1][i-1] += b * Gerror;
                Berr[j+1][i-1] += b * Berror;
            }
            if(j < ht-1){
                Rerr[j+1][i] += x * Rerror;
                Gerr[j+1][i] += x * Gerror;
                Berr[j+1][i] += x * Berror;
            }
            if(i < wd-1 && j < ht-1){
                Rerr[j+1][i+1] += y * Rerror;
                Gerr[j+1][i+1] += y * Gerror;
                Berr[j+1][i+1] += y * Berror;
            }
            if(i < wd-1){
                Rerr[j][i+1] += a * Rerror;
                Gerr[j][i+1] += a * Gerror;
                Berr[j][i+1] += a * Berror;
            }
            quan->setPixel(i, j, 0, Rcopy);
            quan->setPixel(i, j, 1, Gcopy);
            quan->setPixel(i, j, 2, Bcopy);
        }
    }
    return quan;

}

/*
* nearest neighbor sample
*/
Pixel ip_resample_nearest(Image* src, double x, double y) {
    int xRound = (int)(x + 0.5);
    int yRound = (int)(y + 0.5);
	Pixel myPixel = src -> getPixel(xRound, yRound);

	return myPixel;
}

/*
* bilinear sample
*/

Pixel ip_resample_bilinear(Image* src, double x, double y) {
    int lowerX = (int)x;
    int lowerY = (int)y;
    int upperX = (int)(x + 1);
    int upperY = (int)(y + 1);
    double upperWeightx = (x - (double)lowerX)/(upperX - lowerX);
    double upperWeighty = (y - (double)lowerY)/(upperY - lowerY);
    double ULR,ULG,ULB,URR,URG,URB,LLR,LLG,LLB,LRR,LRG,LRB;
    ULR = src -> getPixel(lowerX, lowerY, 0);
    ULG = src -> getPixel(lowerX, lowerY, 1);
    ULB = src -> getPixel(lowerX, lowerY, 2);
    URR = src -> getPixel(upperX, lowerY, 0);
    URG = src -> getPixel(upperX, lowerY, 1);
    URB = src -> getPixel(upperX, lowerY, 2);
    LLR = src -> getPixel(lowerX, upperY, 0);
    LLG = src -> getPixel(lowerX, upperY, 1);
    LLB = src -> getPixel(lowerX, upperY, 2);
    LRR = src -> getPixel(upperX, upperY, 0);
    LRG = src -> getPixel(upperX, upperY, 1);
    LRB = src -> getPixel(upperX, upperY, 2);
    double RmiddleU, GmiddleU, BmiddleU, RmiddleB, GmiddleB, BmiddleB;
    RmiddleU = ULR*(1 - upperWeightx) + URR * upperWeightx;
    GmiddleU = ULG*(1 - upperWeightx) + URG * upperWeightx;
    BmiddleU = ULB*(1 - upperWeightx) + URB * upperWeightx;
    RmiddleB = LLR*(1 - upperWeightx) + LRR * upperWeightx;
    GmiddleB = LLG*(1 - upperWeightx) + LRG * upperWeightx;
    BmiddleB = LLB*(1 - upperWeightx) + LRB * upperWeightx;
    double finalR, finalG, finalB;
    finalR = RmiddleU * (1 - upperWeighty) + RmiddleB * upperWeighty;
    finalG = GmiddleU * (1 - upperWeighty) + GmiddleB * upperWeighty;
    finalB = BmiddleU * (1 - upperWeighty) + BmiddleB * upperWeighty;
    
	Pixel myPixel(finalR,finalG,finalB);
	return myPixel;
}

/*
* gausian sample
*/
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma)
{
    int wd = src -> getWidth();
    int ht = src -> getHeight();
    int startX = (int)x - (size/2 - 1);
    int startY = (int)y - (size/2 - 1);
    double srcR, srcG, srcB, sumR=0, sumG=0, sumB=0, sumWeight=0, weight;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            if(startX + i >= 0 && startX + i < wd && startY + j >=0 && startY + j < ht){
                srcR = src -> getPixel(startX + i, startY + j, 0);
                srcG = src -> getPixel(startX + i, startY + j, 1);
                srcB = src -> getPixel(startX + i, startY + j, 2);
                weight = (1.0/pow(M_E,(sqr(x - (startX + i))+sqr(y - (startY + j)))/(2*sqr(sigma)))) ;
                sumR += weight * srcR;
                sumG += weight * srcG;
                sumB += weight * srcB;
                sumWeight += weight;
            }
        }
    }
    double finalR, finalG, finalB;
    finalR = sumR / sumWeight;
    finalG = sumG / sumWeight;
    finalB = sumB / sumWeight;
	Pixel myPixel(finalR,finalG,finalB);
	return myPixel;
}

/*
* rotate image using one of three sampling techniques
*/
Image* ip_rotate (Image* src, double theta, int x, int y, int samplingMode, 
				  int gaussianFilterSize, double gaussianSigma)
{
    int wd = src -> getWidth();
    int ht = src -> getHeight();
    int bt = src -> getBits();
    double rad = M_PI*2 - deg2rad(theta);
    Image *rotate = new Image(wd,ht,bt);
    double offsetX, offsetY, newX, newY, sampleX,sampleY;
    Pixel pix;
    for(int i = 0 ; i < wd; i++){
        for (int j = 0 ; j < ht; j++){
            offsetX = i - x;
            offsetY = j - y;
            newX = offsetX*cos(rad) - offsetY*sin(rad);
            newY = offsetX*sin(rad) + offsetY*cos(rad);
            sampleX = x + newX;
            sampleY = y + newY;
            if(sampleX < 0 || sampleX > wd-1 || sampleY < 0 || sampleY > ht-1){
                pix = Pixel(0,0,0);
            }else{
                if (samplingMode == 0){
                    pix = ip_resample_nearest(src, sampleX, sampleY);
                }else if(samplingMode == 1){
                    pix = ip_resample_bilinear(src, sampleX, sampleY);
                }else{
                    pix = ip_resample_gaussian(src, sampleX, sampleY, gaussianFilterSize, gaussianSigma);
                }
            }
            rotate -> setPixel(i, j, pix);
        }
    }
	return rotate;
}


/*
* change saturation
*/
Image* ip_saturate (Image* src, double alpha)
{
    Image* greyScale = ip_grey(src);
    return ip_interpolate(src, greyScale, alpha);
}


/*
* scale image using one of three sampling techniques
*/
Image* ip_scale (Image* src, double xFac, double yFac, int samplingMode, 
				 int gaussianFilterSize, double gaussianSigma)
{
    int origWd = src -> getWidth();
    int origHt = src -> getHeight();
    int wd = (src -> getWidth()) * xFac;
    int ht = (src -> getHeight()) * yFac;
    int bt = src -> getBits();
    double originalX, originalY;
    Pixel pix;
    Image *scale = new Image(wd,ht,bt);
    for(int i = 0; i < wd; i++){
        for(int j = 0; j < ht; j++){
            originalX = i / xFac;
            originalY = j / yFac;
            if(originalX < 0 || originalX > origWd-1 || originalY < 0 || originalY > origHt-1){
                pix = Pixel(0,0,0);
            }else{
                if(samplingMode == 0){
                    pix = ip_resample_nearest(src, originalX, originalY);
                }else if(samplingMode == 1){
                    pix = ip_resample_bilinear(src, originalX, originalY);
                }else{
                    pix = ip_resample_gaussian(src, originalX, originalY, gaussianFilterSize, gaussianSigma);
                }
            }
            scale -> setPixel(i, j, pix);
        }
    }
	return scale;
}


/*
* threshold image
*/
Image* ip_threshold (Image* src, double cutoff)
{
    int wd = src->getWidth();
    int ht = src->getHeight();
    int bt = src->getBits();
    Image* threshold = new Image(wd,ht,bt);
    double srcR,srcG,srcB;
    for(int i=0;i<wd;i++){
        for(int j=0;j<ht;j++){
            srcR = src->getPixel(i, j, 0);
            srcG = src->getPixel(i, j, 1);
            srcB = src->getPixel(i, j, 2);
            if(srcR > cutoff){
                threshold->setPixel(i, j, 0, 1);
            }else{
                threshold->setPixel(i, j, 0, 0);
            }
            if(srcG > cutoff){
                threshold->setPixel(i, j, 1, 1);
            }else{
                threshold->setPixel(i, j, 1, 0);
            }
            if(srcB > cutoff){
                threshold->setPixel(i, j, 2, 1);
            }else{
                threshold->setPixel(i, j, 2, 0);
            }
        }
    }
    return threshold;
}




