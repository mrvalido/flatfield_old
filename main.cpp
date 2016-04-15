/*This is the function that does all the dirty work for you, except
;  calculating the displacement between images. This can be done with
;  the routine get_disp included in this package, or use your own
;  favorite method.
;
;  The required input for this function are:
;
;		data:  3-D array containing the succesive displaced images. JAB

;
;		disp:		a 2 x no_of_image real array containning the
;					displacement of each image. The relative
;					displacement will be calculated as
;					dr(i->j) = disp(*,j)-disp(*,i)
					(dx,dy)
;
;		rmin & rmax:Intensity threshold. If a pixel has intensity
;					greater than rmax or less than rmin,
;				        then that pixel will be marked 'bad'
;                                       and will not be used for the gain
;				        calculation.
;
;  On exit, the function returns the following variables:
;
;		flat:		The flatfield after no_iter iterations.
;
;		con:		The algorithm constant (see the paper by Kuhn,
;			        Lin, & Loranz, 1991 Publication of the
;				Astronomical Society of the Pacific, 103,1097)
;
;		gain:		Log10(flat).
;
;		tmp:		template containning the bad pixel maps for the
;				images
;
;		pix_cnt:	No of pair count. N(x) in Kuhn et al paper.
;
;  We keep the output (con,gain,tmp,pix_cnt) so that we can continue
;  iteration using function 'iterate' without going through get_images
;  and get_con again.
;
;  HISTORY:                                                                JAB
;	    Este programa ha sido esencialmente trabajado por              JAB
;	    Phil Wilberg adaptando el programa fortran original            JAB
;	    hecho por Haosheng Lin.                                        JAB
;	    Modified by Jose A. Bonet on 29 August, 1995                   JAB



*/



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


#define IMIN 0 //.8
#define IMAX 82000 //.0
#define LOOPS 10


#include <opencv2/opencv.hpp>
#include "core/flatfield.hpp"
using namespace cv;



void pinta2(ImageValChar& val,int Dy,int Dx, int indice);

// The library is enclosed in a namespace.

int main(){
//(y,x) ojo con los datos de David Orozco son (x,y)
	const int disp[9][2]={{ 0,0},{-30,294},{180,230},{310, 34},{210,-214},{-8,-340},{-220,-220},{-290,10},{-170,270}};





	string nombreImagen;
	char imageName[] = "./im/im0X.fits";
	vector <ImageValInt> datacube;



	ImageValShort tmp (dimX*dimY);//initial 16 bit mask
	// Leer imagenes desde fichero, guardandolas en el vector de datos

	for(unsigned int i = 0; i < no_of_image; i++) {

		imageName[8] = 48 + i;
		nombreImagen = imageName;

		datacube.push_back(readImageFit(nombreImagen));
		getImages(datacube[i], tmp, IMIN, IMAX, i);
	}



	ImageValDouble pixCnt(0.0, datacube[0].size()); //K&Lin Pixel Count
	ImageValDouble con(0.0, datacube[0].size());//K&Lin Constant

	con = getConst(datacube, tmp, pixCnt, disp);

	ImageValDouble pixCntAux(0.0,dimX*dimY);
	pixCntAux= Max(pixCnt, 1.0);
	ImageValDouble gain(0.0,dimX*dimY);
	gain= con / pixCntAux; //gain is normalized K&L constant

	pixCntAux= Min(pixCnt, 1.0);

#ifdef DEBUG
	cout << "GAIN MAX VVVVVVy min: "  <<  gain.max() << "        " << gain.min() << endl;
	ImageValChar pix=escalado8(gain);
	pinta2(pix,dimX, dimY,4);
	waitKey(0);
#endif

	ImageValDouble flat = iterate(con, gain, tmp, pixCnt,disp, LOOPS);


#ifdef DEBUG
	cout << "CON  MAX VVVVVVy min: "  <<  con.max() << "        " << con.min() << endl;
	cout << "GAIN MAX VVVVVVy min: "  <<  gain.max() << "        " << gain.min() << endl;
	cout << "FLAT  MAX VVVVVVy min: "  <<  flat.max() << "        " << flat.min() << endl;
	pix=escalado8(flat);
	pinta2(pix,dimX, dimY,2);
	waitKey(0);
#endif
	//	//Calculo de la ganancia unitaria
	//	ImageValDouble pixCntAux = Max(pixCnt, 1.0);
	//	ImageValDouble gain = con / pixCntAux;

	// Calculo del flatfield


	//	flat = to16U(flat);
	return 0;
}

//pinta 2 function for debug purpose
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


