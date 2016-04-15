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

#define no_of_image 9
#define dimX					2048
#define dimY					2048
#define ind( y, x ) ( y*dimX+x )

using namespace cv;
using namespace std;
using namespace CCfits;


ImageValInt readImageFit(string nombreImagen);
void Mask(ImageValInt& data, ImageValShort& tmp,  const int iMin, const int iMax, int index);

template <typename TT>
void pinta(valarray<TT>& val,int Dy,int Dx, int indice);

ImageValChar escalado8(const ImageValDouble& val);
ImageValChar escalado8(const ImageValShort& val);

int* desplazamientos(const int centros[9][2], int imagenQ, int imagenR);

template <typename T>
ImageValDouble ROI(const valarray<T>& val, int dx, int dy);

template <typename T>
void sumROI(valarray<T>& val, const valarray<T>& ROI, int dx, int dy);

unsigned char toUchart( double n );

ImageValDouble getConst(vector<ImageValInt>& data,
		const ImageValShort& tmp,\
		ImageValDouble& pixCnt,\
		const int disp[9][2]);

void doIteration(const ImageValDouble& con,\
		ImageValDouble& gain,\
		const ImageValShort& tmp,\
		const ImageValDouble& pixCnt,\
		const int disp[9][2]);

ImageValDouble iterate(const ImageValDouble& con, \
		ImageValDouble& gain, \
            const ImageValShort& tmp, \
            const ImageValDouble& pixCnt, \
            const int disp[9][2], \
			const unsigned int loops);


#endif
