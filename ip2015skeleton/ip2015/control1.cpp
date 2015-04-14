#include "control.h"
#include "ip.h"
#include "main.h"
#include <stdlib.h>


/*
* IMPORTANT - DO NOT CHANGE THIS FILE - IMPORTANT
*/


enum {
	M_QUIT = 0,
	M_HELP = 1,

	M_FILE_OPEN =2,
	M_FILE_SAVE =3,
	M_FILE_INFO =4,
	M_FILE_REVERT =5,


	M_PROCESS_BLUR_BOX =6,
	M_PROCESS_BLUR_GAUSSIAN =7,
	M_PROCESS_BLUR_TRIANGLE=8,
	M_PROCESS_BRIGHTEN=9,
	M_PROCESS_COLOR_SHIFT=10,
	M_PROCESS_CONTRAST=11,
	M_PROCESS_COMPOSITE=12,
	M_PROCESS_CROP=13,
	M_PROCESS_EDGE_DETECT=14,
	M_PROCESS_EXTRACT=15,
	M_PROCESS_FUN_WARP=16,
	M_PROCESS_GREY=17,
	M_PROCESS_IMAGE_SHIFT=18,
	M_PROCESS_INVERT=19,
	M_PROCESS_MISC=20,
	M_PROCESS_QUANTIZE_SIMPLE=21,
	M_PROCESS_QUANTIZE_ORDERED=22,
	M_PROCESS_QUANTIZE_FLOYD_STEINBERG=23,
	M_PROCESS_ROTATE=24,
	M_PROCESS_SATURATE=25,
	M_PROCESS_SCALE=26,

	M_PROCESS_SET_SAMPLING_BILINEAR=27,
	M_PROCESS_SET_SAMPLING_NEAREST=28,
	M_PROCESS_SET_SAMPLING_GAUSSIAN=29,
	M_PROCESS_THRESHOLD=30,

	M_VIEW_PIXEL_VALUE=31,

	M_LAST_ENUM
} MENU_ITEMS;


int make_menu ()
{
	int file = glutCreateMenu(menu_func);
	glutAddMenuEntry( "Open...",		M_FILE_OPEN);
	glutAddMenuEntry( "Save...",		M_FILE_SAVE);
	glutAddMenuEntry( "Get Image Info",		M_FILE_INFO);
	glutAddMenuEntry( "Revert",		M_FILE_REVERT);

	int blur = glutCreateMenu(menu_func);
	glutAddMenuEntry( "Box...",		M_PROCESS_BLUR_BOX);
	glutAddMenuEntry( "Gaussian...",	M_PROCESS_BLUR_GAUSSIAN);
	glutAddMenuEntry( "Triangle...",	M_PROCESS_BLUR_TRIANGLE);

	int quantize = glutCreateMenu(menu_func);
	glutAddMenuEntry( "Simple...",	M_PROCESS_QUANTIZE_SIMPLE);
	glutAddMenuEntry( "Ordered...",	M_PROCESS_QUANTIZE_ORDERED);
	glutAddMenuEntry( "Floyd-Steinberg...", M_PROCESS_QUANTIZE_FLOYD_STEINBERG);

	int samplingMode = glutCreateMenu(menu_func);
	glutAddMenuEntry("Bilinear sampling", M_PROCESS_SET_SAMPLING_BILINEAR);
	glutAddMenuEntry("Nearest neighbor sampling", M_PROCESS_SET_SAMPLING_NEAREST);
	glutAddMenuEntry("Gaussian resampling...", M_PROCESS_SET_SAMPLING_GAUSSIAN);

	int warp = glutCreateMenu(menu_func);
	glutAddMenuEntry( "Fun...",		M_PROCESS_FUN_WARP);
	glutAddMenuEntry( "Image shift...",	M_PROCESS_IMAGE_SHIFT);
	glutAddMenuEntry( "Rotate...",	M_PROCESS_ROTATE);
	glutAddMenuEntry( "Scale...",	M_PROCESS_SCALE);
	glutAddSubMenu( "Set Sampling Mode", samplingMode);

	int process = glutCreateMenu(menu_func);
	glutAddSubMenu(   "Blur",		blur);
	glutAddMenuEntry( "Brighten...",	M_PROCESS_BRIGHTEN);
	glutAddMenuEntry("Color shift", M_PROCESS_COLOR_SHIFT);
	glutAddMenuEntry( "Composite...",	M_PROCESS_COMPOSITE);
	glutAddMenuEntry( "Contrast...",	M_PROCESS_CONTRAST);
	glutAddMenuEntry( "Crop...",		M_PROCESS_CROP);
	glutAddMenuEntry( "Edge Detect",	M_PROCESS_EDGE_DETECT);
	glutAddMenuEntry( "Extract...",	M_PROCESS_EXTRACT);
	glutAddMenuEntry( "Grey",		M_PROCESS_GREY);
	glutAddMenuEntry( "Invert",		M_PROCESS_INVERT);
	glutAddMenuEntry( "Misc. effects...", M_PROCESS_MISC);
	glutAddSubMenu(   "Quantize",		quantize);
	glutAddMenuEntry( "Saturate...",	M_PROCESS_SATURATE);
	glutAddMenuEntry( "Threshold...",	M_PROCESS_THRESHOLD);
	glutAddSubMenu( "Warp", warp);





	int main = glutCreateMenu(menu_func);
	glutAddSubMenu(   "File",		file);
	glutAddSubMenu(   "Process",		process);
	glutAddMenuEntry( "View pixel value...", M_VIEW_PIXEL_VALUE);
	glutAddMenuEntry( "Help",		M_HELP);
	glutAddMenuEntry( "Quit",		M_QUIT);

	glutAttachMenu(GLUT_RIGHT_BUTTON);

	return main;
}


static inline void checkStream (const istream& in)
{
	if (in.fail())
	{
		cerr << "Fatal error: stream failed!" << endl;
		exit(-1);
	}
}


void menu_func (int value)
{
	// variables used in the switch statement
	char filename[MAX_LINE];

	switch (value)
	{
	case M_QUIT:  // enum #0
		exit(0);
		break;



	case M_HELP:  // enum #1
		menu_help();
		break;



	case M_FILE_OPEN:   // enum #2
		if (!quietMode)
			cerr << "Open file (string - no spaces) : ";
		cin  >> filename;
		checkStream(cin);
		image_load(filename);
		break;


	case M_FILE_SAVE:   // enum #3
		if (!quietMode)
			cerr << "Save as (string - no spaces) : ";
		cin  >> filename;
		checkStream(cin);
		image_save(filename);
		break;


	case M_FILE_INFO:  // enum #4
		image_print_info();
		break;


	case M_FILE_REVERT:  // enum #5
		image_revert();
		break;

	case M_VIEW_PIXEL_VALUE: // enum #31
		{
			if (!currentImage) 
			{
				cerr << "Sorry, no image is loaded." << endl;
				break;
			}
			if (!quietMode)
			{
				cerr << "Current image width and height: " << currentImage->getWidth()-1 << " " 
					<< currentImage->getHeight()-1 << endl;
			}
			int x=getInt("x value of pixel to view");
			int y=getInt("y value of pixel to view");
			if (x<0 || x>=currentImage->getWidth() || y<0 || y>=currentImage->getHeight())
			{
				cerr << "Invalid pixel location." << endl;
				break;
			}
			cerr << "R: " << currentImage->getPixel(x,y,RED);
			cerr << ", G: " << currentImage->getPixel(x,y,GREEN);
			cerr << ", B: " << currentImage->getPixel(x,y,BLUE) << endl;
			break;
		}

	default:
		process_func(value);
	}
	return;
}

void process_func (int value)
{

	Image* resultImage = NULL;
	static int samplingMode = I_NEAREST;
	static int gaussianFilterSize = 3;
	static double gaussianSigma = 1.0;

	//  check if we have an image to process 
	if (!currentImage)
	{
		cerr << "Sorry, no image is loaded!" << endl;
		return;
	}

	switch (value)
	{
	case M_PROCESS_BLUR_BOX:  // enum #6
		{
			int filterSize = getFilterSize();
			if (filterSize>0)
				resultImage = ip_blur_box(currentImage, filterSize);
			break;
		}


	case M_PROCESS_BLUR_GAUSSIAN:  // enum #7
		{
			int filterSize = getFilterSize();
			if (filterSize<=0)
				break;
			double sigma = getPositiveDouble("sigma");
			if (sigma>0)
				resultImage = ip_blur_gaussian(currentImage, filterSize, sigma);
			break;
		}


	case M_PROCESS_BLUR_TRIANGLE:  // enum #8
		{
			int filterSize = getFilterSize();
			if (filterSize>0)
				resultImage = ip_blur_triangle(currentImage, filterSize);
			break;
		}


	case M_PROCESS_BRIGHTEN:  // enum #9
		{
			double alpha = getDouble("alpha");
			resultImage = ip_brighten(currentImage, alpha);
			break;
		}

	case M_PROCESS_COLOR_SHIFT: // enum #10
		resultImage=ip_color_shift(currentImage);
		break;

	case M_PROCESS_CONTRAST:  // enum #11
		{
			double alpha=getDouble("alpha");
			resultImage = ip_contrast(currentImage, alpha);
			break;
		}


	case M_PROCESS_COMPOSITE: // enum #12
		{
			char filename[MAX_NAME];
			// we don't do a lot of checks here; i.e. second image and
			// mask valid images and the same size as current image
			if (!quietMode)
				cerr << "Enter filename of second image (string - no spaces) : ";
			cin  >> filename;
			Image* secondImage = new Image();
			secondImage->read(filename);
			if (!quietMode)
				cerr << "Enter filename of mask (string - no spaces) : ";
			cin  >> filename;
			Image* mask = new Image();
			mask->read(filename);
			checkStream(cin);
			resultImage = ip_composite(currentImage, secondImage, mask);
			delete secondImage;
			delete mask;
			break;
		}


	case M_PROCESS_CROP: // enum #13
		{
			int x0, y0, x1, y1;
			if (!quietMode)
			{
				cerr << "Current image width and height: " << currentImage->getWidth()-1 << " " 
					<< currentImage->getHeight()-1 << endl;
				cerr << "Enter region to crop (left top right bottom (ints)) : ";
			}
			cin >> x0 >> y0 >> x1 >> y1;
			checkStream(cin);
			if (x0>=0 || x0<x1 || x1<currentImage->getWidth() || y0>=0 || y0<y1 || y1<currentImage->getHeight())
			{
				resultImage = ip_crop(currentImage, x0,y0,x1,y1);
			}
			else
			{
				cerr<< "Invalid region." << endl;
			}

			break;
		}


	case M_PROCESS_EDGE_DETECT: // enum #14
		resultImage = ip_edge_detect(currentImage);
		break;


	case M_PROCESS_EXTRACT:  // enum #15
		{
			int channel = getInt("channel [0,2]");
			if (channel<0 || channel>2) 
			{
				cerr << "Invalid channel."<< endl;
			}
			else
			{
				resultImage = ip_extract(currentImage, channel);
			}
			break;
		}

	case M_PROCESS_FUN_WARP:  // enum #16
		resultImage = ip_fun_warp(currentImage);
		break;

	case M_PROCESS_GREY: // enum #17
		resultImage = ip_grey(currentImage);
		break;

	case M_PROCESS_IMAGE_SHIFT: // enum #18
		{
			double dx = getDouble("dx");
			double dy = getDouble("dy");
			resultImage = ip_image_shift(currentImage,dx, dy);
			break;
		}

	case M_PROCESS_INVERT:  // enum #19
		resultImage = ip_invert(currentImage);
		break;

	case M_PROCESS_MISC: // enum #20
		resultImage = ip_misc(currentImage);
		break;


	case M_PROCESS_QUANTIZE_SIMPLE:  // enum #21
		{

			int bitsPerChannel = getInt("bits per channel [1,8]");
			if (bitsPerChannel<=0 || bitsPerChannel>8)
			{
				cerr << "Invalid number bits." << endl;
			}
			else
			{
				resultImage = ip_quantize_simple(currentImage, bitsPerChannel);
			}
			break;
		}


	case M_PROCESS_QUANTIZE_ORDERED:  //enum #22
		{
			int bitsPerChannel = getInt("bits per channel [1,8]");

			if (bitsPerChannel<=0 || bitsPerChannel>8)
			{
				cerr << "Invalid number bits." << endl;
			}
			else
			{
				resultImage = ip_quantize_ordered(currentImage, bitsPerChannel);
			}
			break;
		}


	case M_PROCESS_QUANTIZE_FLOYD_STEINBERG: // enum #23
		{
			int bitsPerChannel = getInt("bits per channel [1,8]");

			if (bitsPerChannel<=0 || bitsPerChannel>8)
			{
				cerr << "Invalid number bits." << endl;
			}
			else
			{
				resultImage = ip_quantize_fs(currentImage, bitsPerChannel);
			}
			break;
		}


	case M_PROCESS_ROTATE: // enum #24
		{
			if (!quietMode) 
			{
				cerr<< "Current image width/height: " << currentImage->getWidth()-1 << " " 
					<< currentImage->getHeight()-1 << endl;
			}
			double theta = getDouble("angle");
			double x = getDouble("point x");
			double y = getDouble("point y");
			resultImage = ip_rotate(currentImage, theta, x, y, samplingMode, 
				gaussianFilterSize, gaussianSigma);
			break;
		}


	case M_PROCESS_SATURATE:  // enum #25
		{
			double alpha = getDouble("alpha");
			resultImage = ip_saturate(currentImage, alpha);
			break;
		}


	case M_PROCESS_SCALE: // enum #26
		{
			double xFactor = getDouble("scale factor for x");
			double yFactor = getDouble("scale factor for y");
			resultImage = ip_scale(currentImage, xFactor, yFactor, samplingMode, gaussianFilterSize, gaussianSigma);
			break;
		}

	case M_PROCESS_SET_SAMPLING_BILINEAR: // enum #27
		samplingMode=I_BILINEAR;
		break;

	case M_PROCESS_SET_SAMPLING_NEAREST:  // enum #28
		samplingMode=I_NEAREST;
		break;

	case M_PROCESS_SET_SAMPLING_GAUSSIAN: // enum #29
		gaussianFilterSize=getInt("filter size (positive, even integer)");
		if (gaussianFilterSize%2!=0 || gaussianFilterSize<=0)
		{
			gaussianFilterSize=3;
			cerr<<"Invalid value, using default size of 3"<<endl;
		}
		gaussianSigma=getDouble("sigma (non-negative)");
		if (gaussianSigma <0) 
		{
			gaussianSigma=1;
			cerr<<"Invalid value, using default sigma of 1"<<endl;
		}
		samplingMode=I_GAUSSIAN;

	case M_PROCESS_THRESHOLD: // enum #23
		{
			double threshold=getDouble("threshold [0,1]");
			resultImage = ip_threshold(currentImage, threshold);
			break;
		}


	default:
		break;
	}

	if (resultImage != NULL)
	{
		delete currentImage;
		currentImage = resultImage;

		if (currentImage->getWidth()  != window_width    ||
			currentImage->getHeight() != window_height)
			reshape(currentImage->getWidth(), currentImage->getHeight());

		if (!quietMode)
			cerr << "done!" << endl;

		if (!textMode)
			glutPostRedisplay();
	}
}

void keyboard_func (unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'H':
	case 'h':
		menu_help();
		break;
		;;

	case 'Q':
	case 'q':
		exit(0);
		break;
		;;
	}
}


void menu_help ()
{
	cerr << endl
		<< "hmc cs155 image processor" << endl
		<< "please see the ip manual for usage and algorithm information" << endl
		<< "http://www.cs.hmc.edu/courses/2002/fall/cs155/proj1/doc/ip_manual.html"
		<< endl << endl;
}


#define MENUOP(num, tag)	cerr << " " << num << ") " << tag << endl;



void textMenuLoop ()
{
	char command[MAX_LINE];


	while (true)
	{
		if (!quietMode)
			cerr << endl
			<< "selection > " << flush;
		cin  >> command;

		switch (command[0])
		{
		case '\n':
		case '\0':
			//printMenu();
			break;

		case 'Q':
		case 'q':
			menu_func(M_QUIT);
			break;

		case 'H':
		case 'h':
			menu_func(M_HELP);
			break;

		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			menu_func(atoi(command));
			break;

		default:
			//printMenu();
			break;
		}
	}
}


void image_load (const char* filename)
{
	if (currentImage)
		delete currentImage;
	if (originalImage)
		delete originalImage;
	currentImage  = NULL;
	originalImage = NULL;

	originalImage = new Image();
	originalImage->read(filename);

	if (originalImage->good())
	{  
		currentImage = new Image(*originalImage);
		reshape(currentImage->getWidth(), currentImage->getHeight());
	}
	else
	{
		delete originalImage;  
		originalImage = NULL;
		cerr << "Couldn't load image " << filename << "!" << endl;
		return;
	}

	if (!textMode)
		glutPostRedisplay();

	if (!quietMode)
		cerr << "done!" << endl;
}  


void image_save (const char* filename)
{
	if (currentImage)
	{
		if (currentImage->write(filename) == 0)
		{
			//delete originalImage;
			//originalImage = new Image(*currentImage);
		}
	}  
	else if (originalImage)
	{
		originalImage->write(filename);
	}
	else
	{
		cerr << "No image!" << endl;
		return;
	}

	if (!quietMode)
		cerr << "done!" << endl;
}


void image_print_info ()
{  
	if (currentImage) {
		cerr << "width:    " << currentImage->getWidth() << endl
			<< "height:   " << currentImage->getHeight() << endl
			<< "bits:     " << currentImage->getBits() << endl;
	}
	cerr << "done!" << endl;
}


void image_revert ()
{
	if (currentImage)
		delete currentImage;

	if (originalImage)
	{
		currentImage = new Image(*originalImage);

		if (window_width  != currentImage->getWidth() ||
			window_height != currentImage->getHeight())
			reshape(currentImage->getWidth(), currentImage->getHeight());
	}
	else
	{
		cerr << "No image!" << endl;
		return;
	}

	if (!textMode)
		glutPostRedisplay();

	if (!quietMode)
		cerr << "done!" << endl;
}  

int getFilterSize()
{
	int filtersize;
	if (!quietMode)
		cerr << "Enter filter size (positive, odd integer) : ";
	cin  >> filtersize;
	if (filtersize % 2 !=1 || filtersize<=0)
	{
		cerr << "Sorry, the filter size must be a positive, odd integer." << endl;
		filtersize=0;
	}
	checkStream(cin);
	return filtersize;
}
double getDouble(const char* message)
{
	double value;
	if (!quietMode)
		cerr << "Enter " << message << "(double): ";
	cin  >> value;
	checkStream(cin);
	return value;
}

double getPositiveDouble(const char* message)
{
	double value;
	if (!quietMode)
		cerr << "Enter positive " << message << "(double): ";
	cin  >> value;
	checkStream(cin);
	return value;
}

int getInt(const char* message)
{

	int value;
	if (!quietMode)
		cerr << "Enter " << message << " (integer): ";
	cin  >> value;
	checkStream(cin);
	return value;
}
