#include "flatfield.hpp"
#include "utility.hpp"

#define PROGRESS

int Alto;
int Ancho;

ImageValInt readImageFit(string nombreImagen){

	std::auto_ptr<FITS> pInfile(new FITS(nombreImagen,Read,true));
	//std::auto_ptr<FITS> pInfile(new FITS("atestfil.fit",Read,true));

	PHDU& image = pInfile->pHDU();

	valarray<int>  contents;

	// read all user-specifed, coordinate, and checksum keys in the image
	image.readAllKeys();
	image.read(contents);

	int size_val = contents.size();
	ImageValInt im(size_val);
	for(int i = 0; i < size_val; i++){
		im[i] = contents[i];
	}
	return im;
}

void getImages(ImageValInt& data, \
              ImageValChar& tmp, \
              const int iMin, \
              const int iMax,
			  int index) {

	//Obtenemos las dimensiones de una imagen cualquiera (Todas deben ser del mismo tama�o)
	int size_data = data.size();

	// Calcular la plantilla de pixeles buenos
	ImageValChar msk (size_data);
	for(int i=0; i < size_data; i++){
		if(data[i] >= iMin && data[i] <= iMax){
			msk[i] = 1;
		}
	}
	// Guardar la mascara del indezx i en la plantilla de pixeles buenos tmp
	tmp = tmp | (msk * (unsigned char)(1 << index));

	// Emplear la mascara para cribar los pixeles malos
	for(int i=0; i < size_data; i++){
		data[i] = data[i] * (unsigned int)msk[i];
	}
}

ImageValChar escalado8(ImageValDouble& val){
	int size_val = val.size();
	ImageValChar temp(size_val);

	double mx ;
	double min=val.min();
	val=val-min;
	mx=val.max();
	val=val/mx;
	//cout << "Maximo: " << mx << "          Minimo: " << min << endl;

	for(int i = 0; i < size_val; i++){
		temp[i] = (unsigned char) ( val[i] * 255.0);
	}

	return temp;
}
////*************************************************************************************
//ImageValChar escalado8(const ImageValInt& val){
//	int size_val = val.size();
//	ImageValChar temp(size_val);
//
//	unsigned int mx = val.max();
//
//	for(int i = 0; i < size_val; i++){
//		temp[i] = (unsigned char) (( (float)(val[i])/(float)mx ) * 255.0);
//	}
//
//	return temp;
//}

////*************************************************************************************
////*************************************************************************************

int* desplazamientos(int centros[8][2], int imagenQ, int imagenR){
	int* desp = new int[2]; //DY, DX Desplazamiento relativo de las dos imagenes Q y R

	cout << "IMagenQ" << imagenQ  << "      ImagenR: "   << endl;

	desp[0] = centros[imagenQ][0] - centros[imagenR][0];		//DX
	desp[1] = centros[imagenQ][1] - centros[imagenR][1];		//DY

	return desp;
}
template <typename T>
ImageValDouble ROI(const valarray<T>& val, int dx, int dy){

	// Calculo de los extremos de las ventanas
	unsigned int jyl = max(0, -dy);// FILAS
	unsigned int jyh = min(0, -dy) + dimY; // FILAS
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
	//cout << " Alto y ancho : "  <<  alto << "    " << ancho << endl;
    Ancho=ancho;

    Alto=alto;

	ImageValDouble ROI(ta);

	//cout << " Alto: "  <<  ancho << "    " << alto << "    " << ancho*alto << "   "<< ROI.size()<< endl;

	// Calcular ventanas de mascara. MskiqROI y mskirROI son del mismo tamaño, aunque
	//estan desplazadas unas con respecto a la otra una distancia relativa.
	for(int y=jyl; y < jyh; y++){
		for(int x=jxl; x < jxh; x++){
			ROI[(y-jyl)*ancho + (x-jxl)] = (double)val[ind(y,x)];
		}
	}

	return ROI;
}

template <typename T>
void sumROI(valarray<T>& val, valarray<T>& ROI, int dx, int dy){

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


ImageValDouble getConst(vector<ImageValInt>& data, const ImageValChar& tmp, ImageValDouble& pixCnt, const int disp[8][2]) {

	vector<ImageValDouble> dat;
	ImageValDouble con(data[0].size());

	ImageValChar masciq(con.size());//creada para debug
	ImageValChar mascir(con.size());//creada para debug

	// Calculo del logaritmo comun (base 10) de la imagen 0
	dat.push_back(log_10(data[0]));
	//dat=log_10(data[0]);
	int pasada = 0;


	for(unsigned int iq = 1; iq < 8; iq++) {

		// Calculo del logaritmo comun (base 10) de la imagen
		dat.push_back(log_10(data[iq]));
		//dat=log_10(data[iq]);

		// Obtencion de la mascara
		ImageValChar mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			ImageValChar mskir = (tmp & (1 << ir)) / (1 << ir);
			//Desplazamientos
			//int*  desp = desplazamientos(centros, iq, ir);
			//cout << "maximo00  " << (int)mskiq.max() << endl;





			int dx = disp[iq][0] - disp[ir][0];
			int dy = disp[iq][1] - disp[ir][1];
//
//						// Calculo de los extremos de las ventanas
//						unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + data[0].rows; // FILAS
//						unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + data[0].cols; // COLUMNAS
//						unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + data[0].rows; // FILAS
//						unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + data[0].cols; // COLUMNAS



			//Calcula las regiones de interes de  las mascaras
			ImageValDouble mskiqROI = ROI(mskiq, dx, dy);//dx y dy en este

//			mskiqROI=mskiqROI*255;
//			pinta(mskiqROI,Alto,Ancho, 1);
//			waitKey(0);
			ImageValDouble mskirROI = ROI(mskir, -dx, -dy);

//			mskirROI=mskirROI*255;
//			pinta(mskirROI,Alto,Ancho, 2);
//			waitKey(0);
			//sumROI(masciq, mskiq, dx, dy);


			//-----------------------------
			//SOLO DEBUG!!
//			ImageValChar valiq(255,(int)(Alto*Ancho));
//			sumROI(masciq, valiq, dx, dy);
//
//			pinta(masciq,dimX,dimY, 1);
//			waitKey(0);


			//------------------------------

			ImageValDouble mskDouble = mskiqROI * mskirROI;




			//Calcula las regiones de interes de  las mascaras
			ImageValDouble datiqROI = ROI(dat[iq], dx, dy);
			ImageValDouble datirROI = ROI(dat[ir], -dx, -dy);


			ImageValDouble diff = (datiqROI - datirROI)*(mskDouble);


//			ImageValDouble conJROI = ROI(con, dx, dy);
//			ImageValDouble conIROI = ROI(con,  -dx, -dy);
			ImageValDouble conJROI (diff.size());
			ImageValDouble conIROI (diff.size());

			// Aplicar la diferencia a las ventanas del termino constante
			conJROI = diff;
			sumROI(con, conJROI, dx, dy);
			conIROI = - diff;
			sumROI(con, conIROI,  -dx, -dy);


			// Calcular ventanas de la matriz de conteo de pares de pixeles
			ImageValDouble pixCntJROI (mskDouble.size());
			ImageValDouble pixCntIROI (mskDouble.size());

			// Aplicar la mascara a las ventanas de la matriz de pares de pixeles
			pixCntJROI =  mskDouble;
			sumROI(pixCnt, pixCntJROI, dx, dy);
			pixCntIROI =  mskDouble;
			sumROI(pixCnt, pixCntIROI,  -dx, -dy);

			pasada++;
			cout << "Pasada: " << pasada << endl;

		}

	}


	return con;

}

//********************************************************

void doIteration(const ImageValDouble& con,\
		ImageValDouble& gain,\
		const ImageValChar& tmp,\
		const ImageValDouble& pixCnt,\
		const int disp[8][2]) {



	//unsigned int loopCnt = 0;

	// Creacion de la ganancia temporal
	ImageValDouble gainTmp=con;

//	con.copyTo(gainTmp);


	//ImageValDouble gainTmp(con.size());

	for(unsigned int iq = 1; iq < 8; iq++) {

		// Obtencion de la mascara
		ImageValChar mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			ImageValChar mskir = (tmp & (1 << ir)) / (1 << ir);

			// Calcula de los desplazamientos relativos
			int dx = disp[iq][0] - disp[ir][0];
			int dy = disp[iq][1] - disp[ir][1];

			//int*  desp = desplazamientos(centros, iq, ir);

			// Calcular ventanas de mascara
			//Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
			//Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));
			ImageValDouble mskiqROI = ROI(mskiq, dx, dy);
			ImageValDouble mskirROI = ROI(mskir, -dx,-dy);

			// Calcular la mascara de las ventanas
			ImageValDouble mskDouble = mskiqROI * mskirROI;
		//	cout << "GGGGGGGGGGGGGG" << mskDouble.min() << "  "<<mskDouble.max() <<endl;
			//waitKey(0);

//			 msk=(mskiq(jxl:jxh,jyl:jyh) and mskir(ixl:ixh,iyl:iyh))

//			      gain_n(jxl:jxh,jyl:jyh) = gain_n(jxl:jxh,jyl:jyh) $
//			                              + gain(ixl:ixh,iyl:iyh)*msk
//			      gain_n(ixl:ixh,iyl:iyh) = gain_n (ixl:ixh,iyl:iyh) $
//			                              + gain(jxl:jxh,jyl:jyh)*msk

//			Mat msk = mskiqROI.mul(mskirROI);
//			msk.convertTo(msk, CV_64F);

			// Calcular ventanas de ganancia y ganancia temporal
			ImageValDouble gainTmpJROI(mskDouble.size());
			ImageValDouble gainTmpIROI(mskDouble.size());

			ImageValDouble gainJROI=ROI(gain, dx, dy);
			ImageValDouble gainIROI=ROI(gain, -dx,-dy);

			// Modificar la ganancia temporal en base a la ganancia y la mascara
			gainTmpJROI = gainJROI*mskDouble;
			sumROI(gainTmp, gainTmpJROI, dx, dy);

			gainTmpIROI = gainIROI*mskDouble;
			sumROI(gainTmp, gainTmpIROI, -dx,-dy);



#ifdef PROGRESS

			cout << "doItera : Iteración " << " de 28..." << endl;

#endif

		}

	}

	// Calcular ganancia unitaria
	ImageValDouble pixCntAux;
	//Mat pixCntAux = max(pixCnt, 1.0);
	//pixCntAux.convertTo(pixCntAux, CV_64F);
	pixCntAux = Max(pixCnt,1);//normalizo gainTmp
	gainTmp = gainTmp / pixCntAux;

	// Eliminar elementos a cero (de la matriz de pares de pixeles)
//	Mat index = min(pixCnt, 1.0);/marco los pixeles que tiene contribucion para calcular la media
//	index.convertTo(index, CV_64F);

    ImageValDouble ind=Min(pixCnt,1.0);

	//gainTmp = gainTmp.mul(index);
    cout << "minimo GainTM  "<< gainTmp.min() <<"   "<< gainTmp.max() <<endl;
    gainTmp=gainTmp*ind;
    cout << "minimo ind  "<< ind.min() <<"   "<< ind.max() <<endl;
    cout << "minimo GainTM2  "<< gainTmp.min() <<"   "<< gainTmp.max() <<endl;
    //ImageValChar pix=escalado8(gainTmp);
      //  cout  << "idex.........." << endl;
        //pinta(pix,dimX,dimY,1);


    ImageValDouble TMP;
    TMP=gainTmp*gainTmp;

	// Calcular sumatorios

    double sum2=gainTmp.sum();
    double sum3=TMP.sum();
    double nPix=ind.sum();
    cout << "npix     "<< nPix <<  endl;
    cout << " total "<< dimX*dimY <<endl;
    //waitKey(0);
	//double sum2 = sum(gainTmp);
	//double sum3 = sum(gainTmp.mul(gainTmp))[0];
	//double nPix = sum(index)[0];

	// Eliminar elementos mas de 5-Sigma veces alejados de la media
	double ave2 = sum2 / nPix;
	cout<< "media   "<< ave2 << endl;

	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);

	ind = (abs(gainTmp - ave2) > fiveSigma) / 255;
	//index.convertTo(index, CV_64F);
    TMP=gainTmp*ind;
    sum2=sum2-TMP.sum();
    nPix=nPix-ind.sum();
	//sum2 = sum2 - sum(gainTmp.mul(index))[0];
	//nPix = nPix - sum(index)[0];

	// Normalizar la tabla de ganancias
	ave2 = sum2 / nPix;

	gainTmp = gainTmp - ave2;

	// Devolver la tabla de ganancias
	gain = gainTmp;


}


ImageValDouble iterate(const ImageValDouble& con, \
		ImageValDouble& gain, \
            const ImageValChar& tmp, \
            const ImageValDouble& pixCnt, \
            const int disp[8][2], \
			const unsigned int loops) {



	for(unsigned int i = 0; i < loops; i++) {

		doIteration(con, gain, tmp, pixCnt,disp);

#ifdef PROGRESS

		cout << "iterate : Iteración " << i + 1 << " de " << loops << "..." << endl;

#endif

	}

	// Calculo de la imagen de flatfield
	//ImageValDouble flat = gain * log(10.0);
	ImageValDouble flat = pow(10.0,gain);
	cout << "P1 " << endl;

	//flat=exp(flat);
    cout << "P2 " << endl;
	ImageValChar tmpAux;
	tmpAux = (tmp > 0) / 255;
	ImageValDouble tmpAuxDouble(tmp.size());
	//std::valarray<bool> comp =(tmp > 0)/255;
//	tmpAux.convertTo(tmpAux, CV_64F);
	tmpAuxDouble=toDouble(tmpAux);
	flat = flat*tmpAuxDouble;
	cout << "P3 " << endl;
	return flat;

}

