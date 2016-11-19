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

	//12     5     4    32    25     9    28    14    18

	//(y,x) ojo con los datos de David Orozco son (x,y)
	const int disp[9][2]={{0,0},{-30,294},{180,230},{310, 34},{210,-214},{-8,-340},{-220,-220},{-290,10},{-170,270}};


	// (y,x)% test_ND08_EVOL_F0.3_NOROT/
	//const int disp[9][2]={{0,0},{0,13},{12,5},{9,-8},{-4,-11},{-12,0},{-5,11},{8,10},{12,-3}};

	// (y,x)% test_ND08_EVOL_F1_NOROT/
	//const int disp[9][2]={{0,0},{1,147},{ 134,60},{106,-96},{-47,-132},{-142,-5},{-63,132},{94,114},{140,-37}};

	// (y,x)% test_ND08_EVOL_F2_NOROT/
	//const int disp[9][2]={ {0,0},{6,590},{537,241},{427,-384},{-190,-531},{-571,-22},{-254,528},{376,456},{560,-151}};

	//rand 9
	//const int disp[9][2]={{0,0},{1,147},{ 134,60},{427,-384},{-47,-132},{-571,-22},{-254,528},{8,10},{12,-3}};
	//rand 16
	//const int disp[16][2]={{0,0},{12,5},{9,-8},{-12,0},{-5,11},{12,-3},{1,147},{-47,-132},{-63,132},{94,114},{140,-37},{537,241},{427,-384},{-190,-531},{-254,528},{560,-151}};

	//HRT images
	//const int disp[9][2] ={{0,0},{1,147},{134,60},{106,-96},{-47,-132},{-142,-5},{-63,132},{94,114},{140,-37}};

	//HRT disp reales
	//const int disp[9][2] ={{0,0},{2,148},{134,60},{107,-96},{-48,-133},{-143,-6},{-64,132},{94,114},{140,-38}};
//    0           0
//         148           2
//          60         134
//         -96         107
//        -133         -48
//          -6        -143
//         132         -64
//         114          94
//         -38         140

	//Solar C displazamenent
		//const int disp[9][2]={{0,0},{-206,-6},{-178,-106},{-1,-218},{78,-195},{204,4},{98 ,-182},{171,102},{1,206}};

	string nombreImagen;
	//char imageName[] = "./im/im0X.fits";
	//four sets of 9 each displacement image
		char imageName[] = "./im/im0X.fits"; //   8  images set with a displacement below 15% of  solar disc radius
		//char imageName[] = "./imF03/im0X.fits";//11 images set with a displacement around 0.3% of  solar disc radius
		//char imageName[] = "./imF1/im0X.fits";//10  images set with a displacement up to 20% of  solar disc radius
		//char imageName[] = "./imF2/im0X.fits";//  10  images set with a displacement up to 40% of  solar disc radius
	//	char imageName[] = "./imrd/im0X.fits";//  10  images set with a displacement randon 9 solar disc radius
		//char imageName[] = "./im16/imXX.fits";//  10  images set with a displacement randon 16 solar disc radius
		//Test HRT
		//char imageName[] = "./hr/im0X.fits"; //   8  images set HRT

		//TEst Solar C
		//char imageName[] = "./Test_SolarC/im0X.fits"; // 9  images 1080x1080 radio 375
														//desplazamiento maximo 300 entrono al centro, centro 540,540

	vector <ImageValInt> datacube;



	ImageValShort tmp (dimX*dimY);//initial 16 bit mask
	// Leer imagenes desde fichero, guardandolas en el vector de datos

	for(unsigned int i = 0; i < no_of_image; i++) {

		imageName[8] = 48 + i;//for "./im/im0X.fits" set
		//imageName[11] = 48 + i;

		//imageName[17] = 48 + i;//for ./Test_SolarC/im0X.fits" set
		nombreImagen = imageName;

//		if (i<= 9) {
//			imageName[9]=48;
//		}
//		else {
//			imageName[9]=49;

//		}
//		imageName[10] = 48 + i%10;
//		nombreImagen = imageName;
		cout << "nombre image 111111"  	 << endl;
		datacube.push_back(readImageFit(nombreImagen));

		cout << "nombre image"  	 << endl;

		Mask(datacube[i], tmp, IMIN, IMAX, i);
		ImageValChar pixxx=escalado8(datacube[i]);
		pinta2(pixxx,dimX,dimY,i);
		waitKey(0);
	}
	//waitKey(0);

	ImageValDouble pixCnt(0.0, datacube[0].size()); //K&Lin Pixel Count
	ImageValDouble con(0.0, datacube[0].size());//K&Lin Constant

	con = getConst(datacube, tmp, pixCnt, disp);
	ImageValDouble gain(0.0,dimX*dimY);
	//ImageValDouble gain2(0.0,dimX*dimY);
	gain=con;
	//ImageValDouble pixCntAux(0.0,dimX*dimY);
	//pixCntAux= Max(pixCnt, 1.0); //minimo uno
	//gain= con / pixCntAux; //gain is normalized K&L constant
	normalicer(gain, pixCnt);
	cout << "PixCnt MAX VVVVVVy min: "  	<<  pixCnt.max() <<    "     " << pixCnt.min() << endl;
//	cout << "PixCntAux MAX VVVVVVy min: "  	<<  pixCntAux.max() << "     " << pixCntAux.min() << endl;
	cout << "GAIN MAX VVVVVVy min: "  		<<  gain.max() <<      "   	 " << gain.min() << endl;
  //  pixCntAux= Min(pixCnt, 1.0); //imagen de ceros y unos maximo uno

	/////7

	//cout << "PixCntAux MAX VVVVVVy min: "  <<  pixCntAux.max() << "        " << pixCntAux.min() << endl;

	//cout << "GAIN MAX VVVVVVy min: "  <<  gain2.max() << "        " << gain2.min() << endl;

	//



	ImageValDouble flat = iterate(con, gain, tmp, pixCnt,disp, LOOPS);


#ifdef DEBUG
	cout << "CON  MAX VVVVVVy min: "  <<  con.max() << "        " << con.min() << endl;
	cout << "GAIN MAX VVVVVVy min: "  <<  gain.max() << "        " << gain.min() << endl;
	cout << "FLAT  MAX VVVVVVy min: "  <<  flat.max() << "        " << flat.min() << endl;
	ImageValChar pixx=escalado8(con);
	pinta2(pixx,dimX,dimY,1);


	pixx=escalado8(gain);
	pinta2(pixx,dimX,dimY,2);

	ImageValChar pix=escalado8(flat);
	pinta2(pix,dimX, dimY,3);

	int t=writeImage(flat,"flat_hrt.fits",DOUBLE_IMG);
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





///*This is the function that does all the dirty work for you, except
//;  calculating the displacement between images. This can be done with
//;  the routine get_disp included in this package, or use your own
//;  favorite method.
//;
//;  The required input for this function are:
//;
//;		data:  3-D array containing the succesive displaced images. JAB
//
//;
//;		disp:		a 2 x no_of_image real array containning the
//;					displacement of each image. The relative
//;					displacement will be calculated as
//;					dr(i->j) = disp(*,j)-disp(*,i)
//					(dx,dy)
//;
//;		rmin & rmax:Intensity threshold. If a pixel has intensity
//;					greater than rmax or less than rmin,
//;				        then that pixel will be marked 'bad'
//;                                       and will not be used for the gain
//;				        calculation.
//;
//;  On exit, the function returns the following variables:
//;
//;		flat:		The flatfield after no_iter iterations.
//;
//;		con:		The algorithm constant (see the paper by Kuhn,
//;			        Lin, & Loranz, 1991 Publication of the
//;				Astronomical Society of the Pacific, 103,1097)
//;
//;		gain:		Log10(flat).
//;
//;		tmp:		template containning the bad pixel maps for the
//;				images
//;
//;		pix_cnt:	No of pair count. N(x) in Kuhn et al paper.
//;
//;  We keep the output (con,gain,tmp,pix_cnt) so that we can continue
//;  iteration using function 'iterate' without going through get_images
//;  and get_con again.
//;
//;  HISTORY:                                                                JAB
//;	    Este programa ha sido esencialmente trabajado por              JAB
//;	    Phil Wilberg adaptando el programa fortran original            JAB
//;	    hecho por Haosheng Lin.                                        JAB
//;	    Modified by Jose A. Bonet on 29 August, 1995                   JAB
//
//
//
//*/
//
//
//
//// cookbook CCfits demonstration program
////	Astrophysics Science Division,
////	NASA/ Goddard Space Flight Center
////	HEASARC
////	http://heasarc.gsfc.nasa.gov
////	e-mail: ccfits@legacy.gsfc.nasa.gov
////
////	Original author: Ben Dorman
//
//
//// The CCfits headers are expected to be installed in a subdirectory of
//// the include path.
//
//// The <CCfits> header file contains all that is necessary to use both the CCfits
//// library and the cfitsio library (for example, it includes fitsio.h) thus making
//// all of cfitsio's symbolic names available.
//
//#ifdef _MSC_VER
//#include "MSconfig.h" // for truncation warning
//#endif
//
//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif
//
//// this includes 12 of the CCfits headers and will support all CCfits operations.
//// the installed location of the library headers is $(ROOT)/include/CCfits
//
//// to use the library either add -I$(ROOT)/include/CCfits or #include <CCfits/CCfits>
//// in the compilation target.
//
//#define DEBUG
//
//
//#define IMIN 0 //.8
//#define IMAX 82000 //.0
//#define LOOPS 10
//
//
//#include <opencv2/opencv.hpp>
//#include "core/flatfield.hpp"
//using namespace cv;
//
//
//
//void pinta2(ImageValChar& val,int Dy,int Dx, int indice);
//
//// The library is enclosed in a namespace.
//
//int main(){
//
//
//
//	//(y,x) ojo con los datos de David Orozco son (x,y)
//	//const int disp[9][2]={{0,0},{-30,294},{180,230},{310, 34},{210,-214},{-8,-340},{-220,-220},{-290,10},{-170,270}};
//
//	// (y,x)% test_ND08_EVOL_F0.3_NOROT/
//	//const int disp[9][2]={{0,0},{0,13},{12,5},{9,-8},{-4,-11},{-12,0},{-5,11},{8,10},{12,-3}};
//
//	// (y,x)% test_ND08_EVOL_F1_NOROT/
//	//const int disp[9][2]={{0,0},{1,147},{ 134,60},{106,-96},{-47,-132},{-142,-5},{-63,132},{94,114},{140,-37}};
//
//	// (y,x)% test_ND08_EVOL_F2_NOROT/
//	//const int disp[9][2]={ {0,0},{6,590},{537,241},{427,-384},{-190,-531},{-571,-22},{-254,528},{376,456},{560,-151}};
//
//	//rand 9
//	//const int disp[9][2]={{0,0},{1,147},{ 134,60},{427,-384},{-47,-132},{-571,-22},{-254,528},{8,10},{12,-3}};
//	//rand 16
//	const int disp[16][2]={{0,0},{12,5},{9,-8},{-12,0},{-5,11},{12,-3},{1,147},{-47,-132},{-63,132},{94,114},{140,-37},{537,241},{427,-384},{-190,-531},{-254,528},{560,-151}};
//
//	/*
//	//const int disp[9][2]={{  0     0},
//     {1   147},
//   {134    60},
//   {106   -96},
//   {-47  -132},
//  {-142    -5},
//  { -63   132},
//   { 94   114},
//   {140   -37}};
//	 *
//	 * */
//
//	string nombreImagen;
//	//char imageName[] = "./im/im0X.fits";
//	//four sets of 9 each displacement image
//		//char imageName[] = "./im/im0X.fits"; //   8  images set with a displacement below 15% of  solar disc radius
//		//char imageName[] = "./imF03/im0X.fits";//11 images set with a displacement around 0.3% of  solar disc radius
//		//char imageName[] = "./imF1/im0X.fits";//10  images set with a displacement up to 20% of  solar disc radius
//		//char imageName[] = "./imF2/im0X.fits";//  10  images set with a displacement up to 40% of  solar disc radius
//		//char imageName[] = "./imrd/im0X.fits";//  10  images set with a displacement randon 9 solar disc radius
//		//char imageName[] = "./im16/imXX.fits";//  10  images set with a displacement randon 16 solar disc radius
//		//Imagenes HRT
//		char imageName[] = "./imHR/im0X.fits"; //   10
//	vector <ImageValInt> datacube;
//
//
//
//	ImageValShort tmp (dimX*dimY);//initial 16 bit mask
//	// Leer imagenes desde fichero, guardandolas en el vector de datos
//
//	for(unsigned int i = 0; i < no_of_image; i++) {
//
//		if (i<= 9) {
//			imageName[9]=48;
//		}
//		else {
//			imageName[9]=49;
//		}
//		imageName[10] = 48 + i%10;
//		nombreImagen = imageName;
//
//		datacube.push_back(readImageFit(nombreImagen));
//
//		Mask(datacube[i], tmp, IMIN, IMAX, i);
//
//	}
//
//
//	ImageValDouble pixCnt(0.0, datacube[0].size()); //K&Lin Pixel Count
//	ImageValDouble con(0.0, datacube[0].size());//K&Lin Constant
//
//	//con = getConst(datacube, tmp, pixCnt, disp);
//	ImageValDouble gain(0.0,dimX*dimY);
//
//	gain=con;
//
//	normalicer(gain, pixCnt);
//	ImageValDouble flat = iterate(con, gain, tmp, pixCnt,disp, LOOPS);
//
//
//#ifdef DEBUG
//	cout << "CON  MAX VVVVVVy min: "  <<  con.max() << "        " << con.min() << endl;
//	cout << "GAIN MAX VVVVVVy min: "  <<  gain.max() << "        " << gain.min() << endl;
//	cout << "FLAT  MAX VVVVVVy min: "  <<  flat.max() << "        " << flat.min() << endl;
//	ImageValChar pixx=escalado8(con);
//	pinta2(pixx,dimX,dimY,1);
//	pixx=escalado8(gain);
//	pinta2(pixx,dimX,dimY,2);
//	ImageValChar pix=escalado8(flat);
//	pinta2(pix,dimX, dimY,3);
//	int t=writeImage(flat,"flat2.fit",DOUBLE_IMG);
//	waitKey(0);
//#endif
//
//	return 0;
//}
//
////pinta 2 function for debug purpose
//void pinta2(ImageValChar& val,int Dy,int Dx, int indice){
//
//	Mat im(Dy, Dx, CV_8U, Scalar(0));  //Es un tipo de dato de 4 bytes 32S
//
//
//	//Se pone primero el eje Y y despues el eje XCV_64F
//	for (int y=0; y<Dy; y++){
//
//		for (int x=0; x<Dx; x++){
//			//cout << " y   x  : "  << y*Dx + x << "  " << x << endl;
//			im.at<uchar>(y,x) = val[y*Dx + x];
//		}
//	}
//	char imageName[] = "imX.jpg";
//	imageName[2] = 48 + indice;
//	imwrite(imageName, im);
//
//	namedWindow("PINTA", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
//	imshow("PINTA", im);
//}
//
//
