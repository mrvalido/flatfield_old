#include "flatfield.hpp"
#include "utility.hpp"

#define PROGRESS

int Alto;
int Ancho;

ImageValChar masciq(dimX*dimY);//creada para debug
ImageValChar mascir(dimX*dimY);//creada para debug





ImageValInt readImageFit(string nombreImagen){

	std::auto_ptr<FITS> pInfile(new FITS(nombreImagen,Read,true));
	//std::auto_ptr<FITS> pInfile(new FITS("atestfil.fit",Read,true));

	PHDU& image = pInfile->pHDU();

	valarray<int>  contents;

	// read all user-specifed, coordinate, and checksum keys in the image
	image.readAllKeys();
	image.read(contents);
    cout<< "datos de imagen"<< image<<endl;
	int size_val = contents.size();
	ImageValInt im(size_val);

	for(int i = 0; i < size_val; i++){
		im[i] = contents[i];
	}
	return im;
}
/*************************************************************/
/**
 * this Function  writes image in fit file
 * BITPIX code values for FITS image types
 * BYTE_IMG      8
 * SHORT_IMG    16
 * LONG_IMG     32 n c types  is unsigned Int 2bytes
 * LONGLONG_IMG 64
 * FLOAT_IMG   -32   in c types  is float   4bytes
 * DOUBLE_IMG  -64   in c types is double	8bytes
 * @param val	input image to be write in filename
 * @param fileName	name of Fits file
 * @param bitPix number of bit for image pixs (see above)
 * return error code
 */
//template <typename T>
int writeImage(ImageValDouble val,string fileName,long bitPix)
{

    // Create a FITS primary array containing a 2-D image
    // declare axis arrays.
    long naxis    =   2;
    long naxes[2] = { dimX, dimY };
    //long bitPix=LONG_IMG;

    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    std::auto_ptr<FITS> pFits(0);

    try
    {
        // overwrite existing file if the file already exists.

        //const std::string fileName("im00_c.fit");

        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.

        pFits.reset( new FITS(fileName , bitPix , naxis , naxes ) );
    }
    catch (FITS::CantCreate)
    {
          // ... or not, as the case may be.
          return -1;
    }

    long nelements(1);
    long  fpixel(1);


    // Find the total size of the array.
    // this is a little fancier than necessary ( It's only
    // calculating naxes[0]*naxes[1]) but it demonstrates  use of the
    // C++ standard library accumulate algorithm.

    nelements = naxes[0]*naxes[1];


    ImageValDouble array(nelements);
    for (int i = 0; i < nelements; ++i)
    {
        array[i] = val[i];
    }

    pFits->pHDU().addKey("BITPIX", bitPix,"bit per pixel");
    pFits->pHDU().addKey("NAXIS",naxis," number of axis ");
    pFits->pHDU().addKey("NAXIS1",naxes[0]," length of axi x ");
    pFits->pHDU().addKey("NAXIS2",naxes[1]," length of axi y ");

    pFits->pHDU().write(fpixel,nelements,array);


    // PHDU's friend ostream operator. Doesn't print the entire array, just the
    // required & user keywords, and is provided largely for testing purposes [see
    // readImage() for an example of how to output the image array to a stream].

    std::cout << pFits->pHDU() << std::endl;

    return 0;
}




//********************************************************
/**
 * Mask fuction, Purpose is to find the valid pixels in the image Starts from dark current and bad pixel. Also, it corrects image Applying a
 * threshold to the image (lower limit < valid pixels < upper limit). Mask should include Active Regions and Faculae(**)
 * @param data  	 data image   Constant algorithm Term
 * @param tmp		initial mask and store masks for each input image
 * @param iMax 		intesity upper limit
 * @param iMin 		intesity lower limit
 * @return data & tmp     mask data and update initial mask
 */
void Mask(ImageValInt& data, \
              ImageValShort& tmp, \
              const int iMin, \
              const int iMax,
			  int index) {

	//Obtenemos las dimensiones de una imagen cualquiera (Todas deben ser del mismo tama�o)
	int size_data = data.size();

	// Calcular la plantilla de pixeles buenos
	ImageValShort msk (size_data);
	for(int i=0; i < size_data; i++){
		if(data[i] > iMin && data[i] <= iMax){
			msk[i] = 1;
		}
	}

	// Guardar la mascara del indezx i en la plantilla de pixeles buenos tmp
	tmp = tmp | (msk * (unsigned short)(1 << index));

	// Emplear la mascara para cribar los pixeles malos
	for(int i=0; i < size_data; i++){
		data[i] = data[i] * (unsigned int)msk[i];

	}
}

ImageValChar escalado8(const ImageValDouble& val){
	int size_val = val.size();
	ImageValChar temp(size_val);
    ImageValDouble tmp=val;
	double mx ;
	double min=tmp.min();
	tmp=tmp-min;
	mx=tmp.max();
	tmp=tmp/mx;
	//cout << "Maximo: " << mx << "          Minimo: " << min << endl;

	for(int i = 0; i < size_val; i++){
		temp[i] = (unsigned char) ( tmp[i] * 255.0);
	}

	return temp;
}
////*************************************************************************************
ImageValChar escalado8(const ImageValShort& val){
	int size_val = val.size();
	ImageValChar temp(size_val);

	unsigned int mx = val.max();

	for(int i = 0; i < size_val; i++){
		temp[i] = (unsigned char) (( (float)(val[i])/(float)mx ) * 255.0);
	}

	return temp;
}

////*************************************************************************************
ImageValChar escalado8(const ImageValInt& val){
	int size_val = val.size();
	ImageValChar temp(size_val);
	ImageValInt val2=val;
	int minimo=val2.min();
		val2=val2-minimo;
		//int maximo=contents.max();
	 int mx = val2.max();

	for(int i = 0; i < size_val; i++){
		temp[i] = (unsigned char) (( (float)(val2[i])/(float)mx ) * 255.0);
	}

	return temp;
}


////*************************************************************************************
////*************************************************************************************

int* desplazamientos(int centros[8][2], int imagenQ, int imagenR){
	int* desp = new int[2]; //DY, DX Desplazamiento relativo de las dos imagenes Q y R

	cout << "IMagenQ" << imagenQ  << "      ImagenR: "   << endl;

	desp[0] = centros[imagenQ][0] - centros[imagenR][0];		//DX
	desp[1] = centros[imagenQ][1] - centros[imagenR][1];		//DY

	return desp;
}

/**
 * Retrieves a Region of Interest (ROI)from input image . First it calculates the lower and upper corner of ROI
 * @param val	input image
 * @param dx	x component of relative displacement
 * @param dy	y component of relative displacement
 * return ROI	calculated ROI
 */
template <typename T>
ImageValDouble ROI(const valarray<T>& val, int dx, int dy){

	// Calculo de los extremos de las ventanas
	unsigned int jyl = max(0, -dy);// FILAS
	unsigned int jyh = min(0,-dy) + dimY; // FILAS
	unsigned int jxl = max(0, -dx); // COLUMNAS
	unsigned int jxh = min(0, -dx) + dimX; // COLUMNAS

//	unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + dimY; // FILAS
//	unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + dimX; // COLUMNAS
//cout << "jyl: " << jyl << "      jxl: " << jxl << "        jyh: " << jyh << "        jxh: " << jxh << endl;

	int ancho = (int)(jxh-jxl);
	int alto  = (int)(jyh-jyl);
	int ta=ancho*alto;
	//Alt= (jxl-jxh);
	//cout << "tamño : "  <<  ta << "    " << ta << endl;
//	cout << " Alto y ancho : "  <<  alto << "    " << ancho << endl;
    Ancho=ancho;

    Alto=alto;

	ImageValDouble ROI(ta);

//	cout << " ancho y alto : "  <<  ancho << "    " << alto << "    " << ancho*alto <<  endl;

	// Calcular ventanas de mascara. MskiqROI y mskirROI son del mismo tamaño, aunque
	//estan desplazadas unas con respecto a la otra una distancia relativa.
	for(int y=jyl; y < jyh; y++){
		for(int x=jxl; x < jxh; x++){
			ROI[(y-jyl)*ancho + (x-jxl)] = (T)val[ind(y,x)];
		}
	}

	return ROI;
}

/**
 * Update a Region of Interest in input image (ROI).
 * First it calculates the lower and upper corner of ROI and add ROI to image val
 * @param val	input image
 * @param dx	x component of relative displacement
 * @param dy	y component of relative displacement
 * return ROI	calculated ROI
 */
template <typename T>
void sumROI(valarray<T>& val, const valarray<T>& ROI, int dx, int dy){

	// Calculo de los extremos de las ventanas
	unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + dimY; // FILAS
	unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + dimX; // COLUMNAS


	int ancho = (int)(jxh-jxl);

	for(int y=jyl; y < jyh; y++){
		for(int x=jxl; x < jxh; x++){
			 val[ind(y,x)] += ROI[(y-jyl)*ancho + (x-jxl)];
		}
	}
}
template <typename TT>
void pinta(valarray<TT>& val,int Dy,int Dx, int indice){
	Mat im(Dy, Dx, CV_8U, Scalar(0));  //Es un tipo de dato de 4 bytes 32S


	//Se pone primero el eje Y y despues el eje XCV_64F
	for (int y=0; y<Dy; y++){

		for (int x=0; x<Dx; x++){
			//cout << " y   x  : "  << y*Dx + x << "  " << x << endl;
			im.at<TT>(y,x) = val[y*Dx + x];
		}
	}
	char imageName[] = "imX.jpg";
	imageName[2] = 48 + indice;
	imwrite(imageName, im);
	namedWindow("PINTA", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
	imshow("PINTA", im);
}





void normalicer(ImageValDouble& GTmp,ImageValDouble& pixCnt2){

	for (int i=0;i<dimX*dimY;i++){
		if (pixCnt2[i]>1)
			GTmp[i]=GTmp[i]/pixCnt2[i];

}
}
void criba_zero(ImageValDouble& GTmp,ImageValDouble& pixCnt2){

		for (int i=0;i<dimX*dimY;i++){
			if (pixCnt2[i]==0){
				GTmp[i]=0;
			}
	}

}

ImageValDouble media(ImageValDouble& GTmp,ImageValDouble& pixCnt2){
	int npix=0;
	double sum2=0.0;
	double sum3=0.0;
	double GT;
	ImageValDouble M(0.0,4);

	for (int i=0;i<dimX*dimY;i++){
		if (pixCnt2[i]>0){
			sum2=sum2+GTmp[i];
			sum3=sum3+GTmp[i]*GTmp[i];
			npix++;
		}
	}
	M[0]=sum2/npix;//return mean
	M[1]=5*sqrt(sum3/npix-M[0]*M[0]);//return five_sigma
	M[2]=npix;//number of pixels
	M[3]=sum2;//
	return M;
}


ImageValDouble  criba_fivesigma(ImageValDouble& val, double aver2, double fiveSigma){
 	 int npix=0;
 	 double sum2=0.0;
 	 ImageValDouble estadistica(2);
 	 for (int i=0;i<dimX*dimY;i++){
 		 if (abs(val[i] - aver2) >fiveSigma){
 			 sum2=sum2+val[i];
 		     npix++;
 		 }
 	 }
 	 //new mean
 	estadistica[0]=sum2;//return mean
 	estadistica[1]=npix;
 	 return estadistica;
  }


/**
 * fiter to
 */
ImageValDouble flatfield(ImageValDouble& val,const ImageValShort& tmp ){
	 //flat=10.^(gain)*(tmp ne 0)
	 ImageValDouble flat(0.0, val.size());
	 for (int i=0;i<dimX*dimY;i++){
		 if (tmp[i] !=0 )
			 flat[i]=pow(10.0,val[i]);
	 }
	 return flat;
 }

 /**
  * get algorithm's constant term it also calculate the valid pixel pairs count
  *
  *
  *nnnnnn
  *
  *
  */
ImageValDouble getConst(vector<ImageValInt>& data, const ImageValShort& tmp, ImageValDouble& pixCnt, const int disp[9][2]) {

	vector<ImageValDouble> dat;
	ImageValDouble con(data[0].size());



	// Calculo del logaritmo comun (base 10) de la imagen 0
	dat.push_back(log_10(data[0]));

	int pasada = 0;
    bool a;
    a=1;
    int b;
    b=(int)a;
    cout << "bool to inter" << b+1 << endl;
	for(unsigned int iq = 1; iq < no_of_image; iq++) {

		// Calculo del logaritmo comun (base 10) de la imagen
		dat.push_back(log_10(data[iq]));
	//	cout << "minimo maximo dat  "<< iq <<"   "<< dat[iq].min() <<"   "<< dat[iq].max() <<endl;
		//dat=log_10(data[iq]);

		// Obtencion de la mascara
		ImageValShort mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			ImageValShort mskir = (tmp & (1 << ir)) / (1 << ir);


			int dy = (disp[iq][0] - disp[ir][0]);
			int dx = (disp[iq][1] - disp[ir][1]);


			//Calcula las regiones de interes de  las mascaras
			ImageValDouble mskiqROI;//(0.0,Alto*Ancho);

			mskiqROI = ROI(mskiq, -dx, -dy);//dx y dy en este


			ImageValDouble mskirROI;//(0.0,Alto*Ancho);

			mskirROI = ROI(mskir,dx, dy);
			ImageValDouble mskDouble = mskiqROI * mskirROI;

			//------------------------------
		//Calcula las regiones de interes de  las mascaras
			ImageValDouble datiqROI = ROI(dat[iq], -dx, -dy);
			ImageValDouble datirROI = ROI(dat[ir], dx, dy);


			ImageValDouble diff = (datiqROI - datirROI)*(mskDouble);
			ImageValDouble conJROI (diff.size());
			ImageValDouble conIROI (diff.size());

			// Aplicar la diferencia a las ventanas del termino constante
			conJROI = diff;
			sumROI(con, conJROI, -dx, -dy);
			conIROI = - diff;
			sumROI(con, conIROI, dx, dy);


			// Calcular ventanas de la matriz de conteo de pares de pixeles
			ImageValDouble pixCntJROI (mskDouble.size());
			ImageValDouble pixCntIROI (mskDouble.size());

			// Aplicar la mascara a las ventanas de la matriz de pares de pixeles
			pixCntJROI =  mskDouble;
			sumROI(pixCnt, pixCntJROI, -dx, -dy);
			pixCntIROI =  mskDouble;
			sumROI(pixCnt, pixCntIROI,  dx, dy);


			pasada++;
			cout << "Pasada: " << pasada << endl;

		}

	}

	return con;

}

//********************************************************
/**
 * Do an iteration
 * @param
 * @param con  	    Constant algorithm Term
 * @param gain 		is log10(flat)
 * @param tmp		store masks for each input image
 * @param  pixCnt 	No of pair count. N(x) in Kuhn et al paper
 * @return gain     updated gain table
 */
void doIteration(const ImageValDouble& con,\
		ImageValDouble& gain,\
		const ImageValShort& tmp,\
		const ImageValDouble& pixCnt,\
		const int disp[9][2]) {



	//unsigned int loopCnt = 0;

	// Creacion de la ganancia temporal
	ImageValDouble gainTmp=con;

//	con.copyTo(gainTmp);


	//ImageValDouble gainTmp(con.size());
int cont=0;
	for(unsigned int iq = 1; iq < no_of_image; iq++) {

		// Obtencion de la mascara
		ImageValShort mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			ImageValShort mskir = (tmp & (1 << ir)) / (1 << ir);

			// Calcula de los desplazamientos relativos
			int dy = disp[iq][0] - disp[ir][0];
			int dx = disp[iq][1] - disp[ir][1];

			//int*  desp = desplazamientos(centros, iq, ir);

			// Calcular ventanas de mascara

			ImageValDouble mskiqROI = ROI(mskiq, -dx, -dy);
			ImageValDouble mskirROI = ROI(mskir, dx, dy);

			// Calcular la mascara de las ventanas
			ImageValDouble mskDouble = mskiqROI * mskirROI;

			// Calcular ventanas de ganancia y ganancia temporal
			ImageValDouble gainTmpJROI(mskDouble.size());
			ImageValDouble gainTmpIROI(mskDouble.size());

			ImageValDouble gainJROI=ROI(gain, -dx, -dy);
			ImageValDouble gainIROI=ROI(gain, dx,dy);

			// Modificar la ganancia temporal en base a la ganancia y la mascara
			gainTmpJROI = gainJROI*mskDouble;
			sumROI(gainTmp, gainTmpJROI, dx, dy);

			gainTmpIROI = gainIROI*mskDouble;
			sumROI(gainTmp, gainTmpIROI, -dx,-dy);


			cont++;
#ifdef PROGRESS

			//cout << "doItera : Iteración " << cont << " de 36..." << endl;

#endif

		}

	}


	// Calcular ganancia unitaria
	ImageValDouble pixCntAux;
	pixCntAux = pixCnt;
	ImageValDouble matrix(0.0,4);
	ImageValDouble matrix2(0.0,2);
	double sum2=0.0;
	int npix=0;
	normalicer(gainTmp,pixCntAux);
	//criba_zero(gainTmp,pixCntAux);
	matrix=media(gainTmp,pixCntAux);
	cout<< "media y 5 Siigna npix y sum2"<< matrix[0] << "  "<< matrix[1] <<"  "<< matrix[2]<<"  "<< matrix[3]<< endl;

	matrix2=criba_fivesigma(gainTmp, matrix[0], matrix[1]);
	cout<< "Sum2 y npix "<< matrix2[0] << "  "<< matrix2[1] << endl;

	//int npix=criba2(gainTmp, matrix[0], matrix[1]);
	sum2=matrix[4]-matrix2[0];
	npix=(int)(matrix[2]-matrix2[1]);
	cout<< "Sum2 y npix despues"<< sum2 << "  "<< npix << endl;



	double aver=sum2/npix;
	cout<< "Nueva media   "<< aver << endl;

	gain = gainTmp-aver;


//	gainTmp = gainTmp/pixCntAux;
//    ImageValDouble indice(0.0,pixCnt.size());// matrix where pixCnt are zero
//	indice=Min(pixCnt,1.0);
//	cout << "minimo GainTM Antes  "<< gainTmp.min() <<"   "<< gainTmp.max() <<endl;
//    cout << "minimo GainTM2  "<< indice.min() <<"   "<< indice.max() <<endl;
//
//
//
//
//
//	 gainTmp=indice*gainTmp;
//	cout << "minimo GainTM Despues "<< gainTmp.min() <<"   "<< gainTmp.max() <<endl;
//	//ImageValChar pixx2=escalado8(gainTmp);
//    //        cout  << "idex.........." << endl;
//   //pinta(pixx2,dimX,dimY,2);
//
//
//  //  cout << "minimo ind  "<< ind.min() <<"   "<< ind.max() <<endl;
//    cout << "minimo GainTM2  "<< indice.min() <<"   "<< indice.max() <<endl;
//
//
//    ImageValDouble TMP;
//    TMP=gainTmp*gainTmp;
//
//	// Calculate mean
//
//    double sum2=gainTmp.sum();
//    double nPix=indice.sum();
//    double ave2 = sum2 / nPix;
//
//    double sum3=TMP.sum();
//
//    //cout << "npix 1    "<< nPix <<  endl;
////    cout << "total "<< dimX*dimY <<endl;
//
//
//	// Eliminar elementos mas de 5-Sigma veces alejados de la media
//
//	ImageValDouble TMP2;
//	TMP2=gainTmp-ave2;
//	TMP2=TMP2*TMP2;
//	double sum4=TMP2.sum();
//	sum4=sum4/nPix;
//	double fiveSigma=5*sqrt(sum4);
//
//	//double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);
//
	//Recalculo el ind

//	nPix=nPix-npix;
//	sum2=gainTmp.sum();
//	ave2=sum2/nPix;

	//G=(abs(gainTmp - ave2) <= fiveSigma);

//	//valarray<double> indice2=indice;
//    //indice=G*indice;
//	//ImageValDouble G=abs(gainTmp - ave2);
//    ImageValChar pixx2=escalado8(indice);
//    cout  << "idex.........." << endl;
//    pinta(pixx2,dimX,dimY,2);
//


	//waitKey(0);
	//index.convertTo(index, CV_64F);
//    TMP=gainTmp*indice;
//    sum2=TMP.sum();
//    nPix=indice.sum();

   //cout << "npix     "<< nPix <<  endl;
    //cout << " ind.sum()   "<< ind.sum() <<endl;
	//sum2 = sum2 - sum(gainTmp.mul(index))[0];
	//nPix = nPix - sum(index)[0];

	// Normalizar la tabla de ganancias
	//ave2 = sum2 / nPix;

	//gainTmp = gainTmp - ave2;

	// Devolver la tabla de ganancias

//	cout<< "media y sma total   "<< ave2 << "  "<< sum2 << endl;





}


ImageValDouble iterate(const ImageValDouble& con, \
		ImageValDouble& gain, \
            const ImageValShort& tmp, \
            const ImageValDouble& pixCnt, \
            const int disp[9][2], \
			const unsigned int loops) {



	for(unsigned int i = 0; i < loops; i++) {

		doIteration(con, gain, tmp, pixCnt,disp);

#ifdef PROGRESS

		cout << "iterate : Iteración " << i + 1 << " de " << loops << "..." << endl;

#endif

	}

	// Calculo de la imagen de flatfield
	//ImageValDouble flat = gain * log(10.0);

	//ImageValDouble flat = pow(10.0,gain);
	ImageValDouble flat=flatfield(gain,tmp );

//
//
//	//flat=exp(flat);
//    cout << "P2 " << endl;
//	ImageValChar tmpAux;
//	tmpAux = (tmp > 0) / 255;
//	ImageValDouble tmpAuxDouble(tmp.size());
//	//std::valarray<bool> comp =(tmp > 0)/255;
////	tmpAux.convertTo(tmpAux, CV_64F);
//	tmpAuxDouble=toDouble(tmpAux);
//	flat = flat*tmpAuxDouble;
//	cout << "P3 " << endl;
	return flat;

}
