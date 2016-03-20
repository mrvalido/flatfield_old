#ifndef FLATFIELD_HPP
#define FLATFIELD_HPP

#include <opencv2/highgui/highgui.hpp>

#include <vector>

#include <math.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstdlib>   // for srand and rand
#include <ctime>     // for time
#include <CCfits/CCfits>
#include <string>
#include <algorithm>
#include <cmath>
#include "utility.hpp"

#define REL8TO64 				72340172838076673
#define dimX					2048
#define dimY					2048
#define ind( y, x ) ( y*dimX+x )

using namespace cv;
using namespace std;
using namespace CCfits;


ImageValInt readImageFit(string nombreImagen);
void getImages(ImageValInt& data, ImageValChar& tmp,  const int iMin, const int iMax, int index);
ImageValChar escalado8(ImageValDouble& val);
ImageValChar escalado8(const ImageValInt& val);
int* desplazamientos(int centros[8][2], int imagenQ, int imagenR);
template <typename T>
ImageValDouble ROI(const valarray<T>& val, int dx, int dy);
template <typename T>
void sumROI(valarray<T>& val, valarray<T>& ROI, int dx, int dy);
ImageValDouble getConst(vector<ImageValInt>& data, const ImageValChar& tmp, ImageValDouble& pixCnt, int centros[8][2]);

void doIteration(const Mat&, Mat&, const Mat&, const Mat&, const int[8][2]);

Mat iterate(const Mat&, Mat&, const Mat&, const Mat&, const int[8][2], const unsigned int);
Mat binarizar (const Mat&);

//NUEVAS FUNCIONES
/*
void getCon(const Mat& mskiq,const Mat& mskir,const Mat& dataiq,const Mat& datair, Mat& pixCnt, Mat& con, int dx, int dy);
Mat getConstNueva(vector<Mat>& data, const Mat& tmp, Mat& pixCnt, int **disp);
void getGainTmp(const Mat& mskiq,const Mat& mskir, Mat& gainTmp, Mat& gain, int dx, int dy);
void calculateStats(const Mat& pixCnt, Mat& gainTmp, Mat& gain);
void doIterationNueva(const Mat& con, Mat& gain, const Mat& tmp, const Mat& pixCnt, const int disp[8][2]);
*/
#endif
