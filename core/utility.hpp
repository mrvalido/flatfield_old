#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <algorithm>
#include <valarray>
using namespace cv;
using namespace std;

typedef valarray<unsigned int>   ImageValInt;
typedef valarray<unsigned char>  ImageValChar;
typedef valarray<unsigned long>  ImageValLong;
typedef valarray<float>  		 ImageValFloat;
typedef valarray<double>		 ImageValDouble;



ImageValDouble log_10(ImageValInt val);

Mat to16U(const Mat&);

#endif
