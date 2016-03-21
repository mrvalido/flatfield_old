#include "flatfield.hpp"
#include "utility.hpp"


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

	double mx = val.max();

	for(int i = 0; i < size_val; i++){
		temp[i] = (unsigned char) (( (double)(val[i])/(double)mx ) * 255.0);
	}

	return temp;
}
////*************************************************************************************
ImageValChar escalado8(const ImageValInt& val){
	int size_val = val.size();
	ImageValChar temp(size_val);

	unsigned int mx = val.max();

	for(int i = 0; i < size_val; i++){
		temp[i] = (unsigned char) (( (float)(val[i])/(float)mx ) * 255.0);
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
	//int alto  = (int)(jyh-jyl);
	//Alt= (jxl-jxh);
	//cout << "Anch : "  <<  ancho << "    " << abs(ancho) << endl;
	//cout << " Alto: "  <<  alto << "    " << abs(alto) << endl;
    //Ancho=ancho;

    //Alto=alto;
	ImageValDouble ROI((jyh-jyl) * (jxh-jxl));

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
			 val[ind(y,x)] = ROI[(y-jyl)*ancho + (x-jxl)];
		}
	}
}

ImageValDouble getConst(vector<ImageValInt>& data, const ImageValChar& tmp, ImageValDouble& pixCnt, const int centros[8][2]) {

	vector<ImageValDouble> dat;

	ImageValDouble con(data[0].size());

	// Calculo del logaritmo comun (base 10) de la imagen
	dat.push_back(log_10(data[0]));

	for(unsigned int iq = 1; iq < 8; iq++) {

		// Calculo del logaritmo comun (base 10) de la imagen
		dat.push_back(log_10(data[iq]));

		// Obtencion de la mascara
		ImageValChar mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			ImageValChar mskir = (tmp & (1 << ir)) / (1 << ir);
			//Desplazamientos
			int*  desp = desplazamientos(centros, iq, ir);

			ImageValDouble mskiqROI = ROI(mskiq, desp[0], desp[1]);
			ImageValDouble mskirROI = ROI(mskir, -desp[0], -desp[1]);

			ImageValDouble mskDouble = mskiqROI * mskirROI;

			ImageValDouble datiqROI = ROI(dat[iq], desp[0], desp[1]);
			ImageValDouble datirROI = ROI(dat[ir], -desp[0], -desp[1]);


			ImageValDouble diff = (datiqROI - datirROI)*(mskDouble);


			ImageValDouble conJROI = ROI(con, desp[0], desp[1]);
			ImageValDouble conIROI = ROI(con, -desp[0], -desp[1]);

			// Aplicar la diferencia a las ventanas del termino constante
			conJROI = conJROI + diff;
			sumROI(con, conJROI, desp[0], desp[1]);
			conIROI = conIROI - diff;
			sumROI(con, conIROI, -desp[0], -desp[1]);


			// Calcular ventanas de la matriz de conteo de pares de pixeles
			ImageValDouble pixCntJROI = ROI(pixCnt, desp[0], desp[1]);
			ImageValDouble pixCntIROI = ROI(pixCnt, -desp[0], -desp[1]);

			// Aplicar la mascara a las ventanas de la matriz de pares de pixeles
			pixCntJROI = pixCntJROI + mskDouble;
			sumROI(pixCnt, pixCntJROI, desp[0], desp[1]);
			pixCntIROI = pixCntIROI + mskDouble;
			sumROI(pixCnt, pixCntIROI, -desp[0], -desp[1]);

		}

	}
	return con;

}

//********************************************************

//void doIteration(const ImageValDouble& con,\
//		ImageValDouble& gain,\
//		const ImageValChar& tmp,\
//		const ImageValDouble& pixCnt,\
//		const int centros[8][2]) {
//
//
//
//	//unsigned int loopCnt = 0;
//
//	// Creacion de la ganancia temporal
////	Mat gainTmp;
////	con.copyTo(gainTmp);
//
//
//	ImageValDouble gainTmp(con.size());
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Obtencion de la mascara
//		ImageValChar mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			ImageValChar mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//		//	int dx = disp[iq][0] - disp[ir][0];
//			//int dy = disp[iq][1] - disp[ir][1];
//			int*  desp = desplazamientos(centros, iq, ir);
////			// Calculo de los extremos de las ventanas
////			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + con.rows; // FILAS
////			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + con.cols; // COLUMNAS
////			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + con.rows; // FILAS
////			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + con.cols; // COLUMNAS
//
//			// Calcular ventanas de mascara
//			//Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
//			//Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));
//			ImageValDouble mskiqROI = ROI(mskiq, desp[0], desp[1]);
//			ImageValDouble mskirROI = ROI(mskir, -desp[0], -desp[1]);
//
//			// Calcular la mascara de las ventanas
//			ImageValDouble mskDouble = mskiqROI * mskirROI;
//
//
//
////			Mat msk = mskiqROI.mul(mskirROI);
////			msk.convertTo(msk, CV_64F);
//
//			// Calcular ventanas de ganancia y ganancia temporal
//			ImageValDouble gainTmpJROI=ROI(gainTmp, desp[0], desp[1]);
//			ImageValDouble gainTmpIROI=ROI(gainTmp, -desp[0], -desp[1]);
//			ImageValDouble gainJROI=ROI(gain, desp[0], desp[1]);
//			ImageValDouble gainIROI=ROI(gain, -desp[0], -desp[1]);
//
//			// Modificar la ganancia temporal en base a la ganancia y la mascara
//			gainTmpJROI = gainTmpJROI + gainIROI*mskDouble;
//			gainTmpIROI = gainTmpIROI + gainJROI*mskDouble;
//
//#ifdef PROGRESS
//
//			cout << "doItera : Iteración " << ++loopCnt << " de 28..." << endl;
//
//#endif
//
//		}
//
//	}
//
////	// Calcular ganancia unitaria
////	ImageValDouble pixCntAux=pixCnt;
////	//Mat pixCntAux = max(pixCnt, 1.0);
////	//pixCntAux.convertTo(pixCntAux, CV_64F);
////	pixCntAux = Max(pixCnt,1);
////	gainTmp = gainTmp / pixCntAux;
////
////	// Eliminar elementos a cero (de la matriz de pares de pixeles)
//////	Mat index = min(pixCnt, 1.0);
//////	index.convertTo(index, CV_64F);
////
////    ImageValDouble index=Min(pixCnt,1);
////	//gainTmp = gainTmp.mul(index);
////    gainTmp=gainTmp*index;
////    ImageValDouble TMP;
////    TMP=gainTmp*gainTmp;
////
////	// Calcular sumatorios
////
////    double sum2=gainTmp.sum();
////    double sum3=TMP.sum();
////    double nPix=index.sum();
////	//double sum2 = sum(gainTmp);
////	//double sum3 = sum(gainTmp.mul(gainTmp))[0];
////	//double nPix = sum(index)[0];
////
////	// Eliminar elementos mas de 5-Sigma veces alejados de la media
////	double ave2 = sum2 / nPix;
////	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);
////
////	index = (abs(gainTmp - ave2) > fiveSigma) / 255;
////	//index.convertTo(index, CV_64F);
////    TMP=gainTmp*index;
////    sum2=sum2-TMP.sum();
////    nPix=nPix-index.sum();
////	//sum2 = sum2 - sum(gainTmp.mul(index))[0];
////	//nPix = nPix - sum(index)[0];
////
////	// Normalizar la tabla de ganancias
////	ave2 = sum2 / nPix;
////
////	gainTmp = gainTmp - ave2;
////
////	// Devolver la tabla de ganancias
////	gain = gainTmp;
//
//}


//ImageValDouble iterate(const ImageValDouble& con, \
//		ImageValDouble& gain, \
//            const ImageValChar& tmp, \
//            const ImageValDouble& pixCnt, \
//            const int centros[8][2], \
//			const unsigned int loops) {
//
//
//
//	for(unsigned int i = 0; i < loops; i++) {
//
//		doIteration(con, gain, tmp, pixCnt,centros);
//
//#ifdef PROGRESS
//
//		cout << "iterate : Iteración " << i + 1 << " de " << loops << "..." << endl;
//
//#endif
//
//	}
//
//	// Calculo de la imagen de flatfield
//	ImageValDouble flat = gain * log(10.0);
//
//	flat=exp(flat);
////
//	ImageValChar tmpAux;
//	//tmpAux = (tmp > 0) / 255;
//	//std::valarray<bool> comp =(tmp > 0)/255;
////	tmpAux.convertTo(tmpAux, CV_64F);
////
//	flat = flat*toDouble(tmpAux);
//
//	return flat;
//
//}

