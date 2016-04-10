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

	int size_val = contents.size();
	ImageValInt im(size_val);

	for(int i = 0; i < size_val; i++){
		im[i] = contents[i];
	}
	return im;
}

void getImages(ImageValInt& data, \
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
	unsigned int jyh = min(0,-dy) + dimY; // FILAS
	unsigned int jxl = max(0, -dx); // COLUMNAS
	unsigned int jxh = min(0, -dx) + dimX; // COLUMNAS

//	unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + dimY; // FILAS
//	unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + dimX; // COLUMNAS

//	cout << "jyl: " << jyl << "      jxl: " << jxl << "        jyh: " << jyh << "        jxh: " << jxh << endl;

	int ancho = (int)(jxh-jxl);
	int alto  = (int)(jyh-jyl);
	int ta=ancho*alto;
	//Alt= (jxl-jxh);
	//cout << "tamño : "  <<  ta << "    " << ta << endl;
	//cout << " Alto y ancho : "  <<  alto << "    " << ancho << endl;
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

 unsigned char toUchart( double n )
{
   return (unsigned char)n;
}

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
	for(unsigned int iq = 1; iq < 9; iq++) {

		// Calculo del logaritmo comun (base 10) de la imagen
		dat.push_back(log_10(data[iq]));
		cout << "minimo maximo dat  "<< iq <<"   "<< dat[iq].min() <<"   "<< dat[iq].max() <<endl;
		//dat=log_10(data[iq]);

		// Obtencion de la mascara
		ImageValShort mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			ImageValShort mskir = (tmp & (1 << ir)) / (1 << ir);
			//Desplazamientos
			//int*  desp = desplazamientos(centros, iq, ir);
			//cout << "maximo00  " << (int)mskiq.max() << endl;





			int dy = (disp[iq][0] - disp[ir][0]);
			int dx = (disp[iq][1] - disp[ir][1]);
//

		//	cout << "dx: " << dx << iq << endl;
		//	cout << "dy: " << dy << ir << endl;


			//Calcula las regiones de interes de  las mascaras
			ImageValDouble mskiqROI(0.0,Alto*Ancho);

			mskiqROI = ROI(mskiq, dx, dy);//dx y dy en este
			ImageValDouble mskirROI(0.0,Alto*Ancho);

			mskirROI = ROI(mskir,-dx, -dy);
			ImageValDouble mskDouble = mskiqROI * mskirROI;
			//-----------------------------
			//SOLO DEBUG!!
//			ImageValChar valiq(50,(Alto*Ancho));
//			cout << "iq   "<< iq<< "ir  "<< ir << endl;
//			masciq=mskiq*50;
//			sumROI(masciq, valiq, -dx, -dy);
//
//			pinta(masciq,dimX,dimY, iq);
//			masciq=mskir*50;
//			sumROI(masciq, valiq, dx, dy);
//			pinta(masciq,dimX,dimY, ir);
//			 waitKey(0);


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
			sumROI(pixCnt, pixCntJROI, dx, dy);
			pixCntIROI =  mskDouble;
			sumROI(pixCnt, pixCntIROI,  -dx, -dy);


			pasada++;
			cout << "Pasada: " << pasada << endl;

		}

	}
	cout << "Pinta COM: " << pasada << endl;
	cout << "COM MAX VVVVVVy min: "  <<  con.max() << "        " <<con.min() << endl;
	masciq=escalado8(con);
	pinta(masciq,dimX,dimY, 2);
	waitKey(0);
	return con;

}

//********************************************************

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

	for(unsigned int iq = 1; iq < 9; iq++) {

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



			//------------------------------

			// Calcular ventanas de ganancia y ganancia temporal
			ImageValDouble gainTmpJROI(mskDouble.size());
			ImageValDouble gainTmpIROI(mskDouble.size());

			ImageValDouble gainJROI=ROI(gain, dx, dy);
			ImageValDouble gainIROI=ROI(gain, -dx,-dy);

			// Modificar la ganancia temporal en base a la ganancia y la mascara
			gainTmpJROI = gainJROI*mskDouble;
			sumROI(gainTmp, gainTmpJROI, -dx, -dy);

			gainTmpIROI = gainIROI*mskDouble;
			sumROI(gainTmp, gainTmpIROI, dx,dy);



#ifdef PROGRESS

			cout << "doItera : Iteración " << " de 28..." << endl;

#endif

		}

	}


	// Calcular ganancia unitaria
	ImageValDouble pixCntAux;
	pixCntAux = Max(pixCnt, 1.0);

	gainTmp = gainTmp/pixCntAux;

    ImageValDouble indice(0.0,pixCnt.size());
	indice=Min(pixCnt,1.0);
	cout << "minimo GainTM Antes  "<< gainTmp.min() <<"   "<< gainTmp.max() <<endl;
    cout << "minimo GainTM2  "<< indice.min() <<"   "<< indice.max() <<endl;

//	ImageValChar pixx=escalado8(gainTmp);
//	pinta(pixx,dimX,dimY,1);
//	waitKey(0);


	//gainTmp = gainTmp.mul(index);


	 gainTmp=indice*gainTmp;
	cout << "minimo GainTM Despues "<< gainTmp.min() <<"   "<< gainTmp.max() <<endl;
//	ImageValChar pixx2=escalado8(gainTmp);
//    //        cout  << "idex.........." << endl;
//   pinta(pixx2,dimX,dimY,2);
//   waitKey(0);

  //  cout << "minimo ind  "<< ind.min() <<"   "<< ind.max() <<endl;
    cout << "minimo GainTM2  "<< indice.min() <<"   "<< indice.max() <<endl;


    ImageValDouble TMP;
    TMP=gainTmp*gainTmp;

	// Calcular sumatorios

    double sum2=gainTmp.sum();
    double sum3=TMP.sum();
    double nPix=indice.sum();

    cout << "npix     "<< nPix <<  endl;
    cout << " total "<< dimX*dimY <<endl;

	//double sum2 = sum(gainTmp);
	//double sum3 = sum(gainTmp.mul(gainTmp))[0];
	//double nPix = sum(index)[0];

	// Eliminar elementos mas de 5-Sigma veces alejados de la media
	double ave2 = sum2 / nPix;
	cout<< "media y sma total   "<< ave2 << "  "<< sum2 << endl;

	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);
	cout<< "cinco sigma   "<< fiveSigma << endl;
	//Recalculo el ind
	valarray<double> G ;
	G=(abs(gainTmp - ave2) <= fiveSigma);
	//valarray<double> indice2=indice;
    indice=G*indice;
	//ImageValDouble G=abs(gainTmp - ave2);
	ImageValChar pixx=escalado8(indice);
	pinta(pixx,dimX,dimY,1);
	waitKey(0);

	cout << "minimo ind  "<< G.min() <<"   "<< G.max() <<endl;

	//waitKey(0);
	//index.convertTo(index, CV_64F);
    TMP=gainTmp*indice;
    sum2=TMP.sum();
    nPix=indice.sum();

  //  cout << "npix     "<< nPix <<  endl;
    //cout << " ind.sum()   "<< ind.sum() <<endl;
	//sum2 = sum2 - sum(gainTmp.mul(index))[0];
	//nPix = nPix - sum(index)[0];

	// Normalizar la tabla de ganancias
	ave2 = sum2 / nPix;

	gainTmp = gainTmp - ave2;

	// Devolver la tabla de ganancias
	gain = gainTmp;
	cout<< "media y sma total   "<< ave2 << "  "<< sum2 << endl;
	pixx=escalado8(gain);
	cout  << "idex.........." <<  endl;
	pinta(pixx,dimX,dimY,3);





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
	ImageValDouble flat = pow(10.0,gain);
	cout << "P1 pinta flat " << endl;
	ImageValChar pixx=escalado8(flat);
	//			//        cout  << "idex.........." << endl;
	pinta(pixx,dimX,dimY,1);
	waitKey(0);
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

