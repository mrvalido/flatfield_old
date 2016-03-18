// cookbook CCfits demonstration program
//	Astrophysics Science Division,
//	NASA/ Goddard Space Flight Center
//	HEASARC
//	http://heasarc.gsfc.nasa.gov
//	e-mail: ccfits@legacy.gsfc.nasa.gov
//
//	Original author: Ben Dorman


// The CCfits headers are expected to be installed in a subdirectory of
// the include path.

// The <CCfits> header file contains all that is necessary to use both the CCfits
// library and the cfitsio library (for example, it includes fitsio.h) thus making
// all of cfitsio's symbolic names available.

#ifdef _MSC_VER
#include "MSconfig.h" // for truncation warning
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// this includes 12 of the CCfits headers and will support all CCfits operations.
// the installed location of the library headers is $(ROOT)/include/CCfits

// to use the library either add -I$(ROOT)/include/CCfits or #include <CCfits/CCfits>
// in the compilation target.



#include <opencv2/opencv.hpp>
#include "core/flatfield.hpp"
using namespace cv;



// The library is enclosed in a namespace.

int main(){

	int centros[8][2] = {
						{1316,991},
				 	    {1252,1201},
						{1056,1331},
						{808,1231},
						{682,1013},
						{802,801},
						{1032,731},
						{1292,851}
					   };


	string nombreImagen;
	char imageName[] = "./im/im0X.fits";
	vector <ImageValInt> datacube;

	int iMin,iMax;

	ImageValChar tmp (dimX*dimY);
	// Leer imagenes desde fichero, guardandolas en el vector de datos
	for(unsigned int i = 0; i < 8; i++) {
		imageName[8] = 49 + i;
		nombreImagen = imageName;
		datacube.push_back(readImageFit(nombreImagen));
		iMin = datacube[i].min();
		iMax = datacube[i].max();

		getImages(datacube[i], tmp, iMin+100, iMax-500, 0);
	}


	ImageValDouble pixCnt(datacube[0].size());

	ImageValDouble con = getConst(datacube, tmp, pixCnt, centros);

	cout << "MAX: "  <<  con.max() << "        " << pixCnt.max() << endl;

	ImageValDouble triow = log_10(datacube[0]);

	ImageValChar im8 = escalado8(triow);

	Mat im(dimY, dimX, CV_8UC1, Scalar(0));  //Es un tipo de dato de 4 bytes 32S

	//Se pone primero el eje Y y despues el eje X
	for (long y=0; y<dimY; y++){
			for (long x=0; x<dimX; x++){
			im.at<uchar>(y,x) = (uchar)(im8[ind( y, x )]);
		}
	}


	imwrite("msk.jpeg", im);
	namedWindow( "Display window", WINDOW_NORMAL);// Create a window for display.
	imshow( "Display window", im );




//	Mat im2(dimY, dimX, CV_8UC1, Scalar(0));  //Es un tipo de dato de 4 bytes 32S
//	//Se pone primero el eje Y y despues el eje X
//	for (long y=0; y<dimY; y++){
//		for (long x=0; x<dimX; x++){
//			im2.at<uchar>(y,x) = (uchar)(open[ind( y, x )]);
//		}
//	}
//
//	namedWindow( "Display window 2", WINDOW_NORMAL );// Create a window for display.
//	imshow( "Display window 2", im2 );

	waitKey(0);
	//imageVal.showImageMat();


	return 0;
}

