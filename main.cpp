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

#define DEBUG


#define IMIN 30096 //.8
#define IMAX 65535 //.0
#define LOOPS 10

#include <opencv2/opencv.hpp>
#include "core/flatfield.hpp"
using namespace cv;



void pinta2(ImageValChar& val,int Dy,int Dx, int indice);

// The library is enclosed in a namespace.

int main(){

//	int centros[8][2] = {
//						{1316,991},
//				 	    {1252,1201},
//						{1056,1331},
//						{808,1231},
//						{682,1013},
//						{802,801},
//						{1032,731},
//						{1292,851}
//					   };

	const int disp[8][2] ={{0, 0},{-64,210},{-260, 340},{-508, 240},{-634,22},{-514,-190},{-284,-260},{-24,-140}};


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
		getImages(datacube[i], tmp, IMIN, IMAX, i);

	}



	ImageValDouble pixCnt(datacube[0].size());
	ImageValDouble con(datacube[0].size());

	con = getConst(datacube, tmp, pixCnt, disp);


	//pinta(con, dimX, dimX, 1);

	cout << "MAX: "  <<  con.max() << "        " << pixCnt.min() << endl;

	ImageValDouble pixCntAux = Max(pixCnt, 1.0);
	cout << "minimo pxaux"<< pixCntAux.min() <<"   "<<pixCntAux.max() <<endl;
	ImageValDouble gain = con / pixCntAux;


	//ImageValChar pix=escalado8(pixCntAux);
	//cout << "minimo pxaux"<< (int) pix.min() <<"   "<< (int)pix.max() <<endl;
	//pinta2(pix,dimX, dimY,1);
	//waitKey(0);
	//cout<< "klsdflsfklskflskfs"<<endl;

	cout << "MAX y min: "  <<  gain.max() << "        " << pixCntAux.max() << endl;
	ImageValDouble flat = iterate(con, gain, tmp, pixCnt,disp, LOOPS);
#ifdef DEBUG


	cout << "GAIN MAX VVVVVVy min: "  <<  gain.max() << "        " << gain.min() << endl;
	cout << "CON  MAX VVVVVVy min: "  <<  con.max() << "        " << con.min() << endl;



	cout << "Pasadaaaaaaaaaaaaaaaaaaa: " << endl;
	//cout << "minimo pxaux"<< (int) pix.min() <<"   "<< (int)pix.max() <<endl;
	ImageValChar pix=escalado8(flat);
	pinta2(pix,dimX, dimY,1);
	waitKey(0);


#endif
	//	//Calculo de la ganancia unitaria
	//	ImageValDouble pixCntAux = Max(pixCnt, 1.0);
	//	ImageValDouble gain = con / pixCntAux;

	// Calculo del flatfield


	//	flat = to16U(flat);
	return 0;
}


void pinta2(ImageValChar& val,int Dy,int Dx, int indice){

	Mat im(Dy, Dx, CV_8U, Scalar(0));  //Es un tipo de dato de 4 bytes 32S


	//Se pone primero el eje Y y despues el eje XCV_64F
	for (int y=0; y<Dy; y++){

		for (int x=0; x<Dx; x++){
			//cout << " y   x  : "  << y*Dx + x << "  " << x << endl;
			im.at<uchar>(y,x) = val[y*Dx + x];
		}
	}
	char imageName[] = "imX.jpg";
	imageName[2] = 48 + indice;
	imwrite(imageName, im);

	namedWindow("PINTA", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
	imshow("PINTA", im);
}


