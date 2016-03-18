//#include "flatfield.hpp"
//#include "utility.hpp"
//
//
//ImageValInt readImageFit(string nombreImagen){
//
//	std::auto_ptr<FITS> pInfile(new FITS(nombreImagen,Read,true));
//	//std::auto_ptr<FITS> pInfile(new FITS("atestfil.fit",Read,true));
//
//	PHDU& image = pInfile->pHDU();
//
//	valarray<int>  contents;
//
//	// read all user-specifed, coordinate, and checksum keys in the image
//	image.readAllKeys();
//	image.read(contents);
//
//	int size_val = contents.size();
//	ImageValInt im(size_val);
//	cout << image << std::endl;
//	for(int i = 0; i < size_val; i++){
//		im[i] = contents[i];
//	}
//	return im;
//}
//
//void getImages(ImageValInt& data, \
//              ImageValChar& tmp, \
//              const int iMin, \
//              const int iMax,
//			  int index) {
//
//	//Obtenemos las dimensiones de una imagen cualquiera (Todas deben ser del mismo tama�o)
//	int size_data = data.size();
//
//	// Calcular la plantilla de pixeles buenos
//	ImageValChar msk (size_data);
//
//	for(int i=0; i < size_data; i++){
//		if(data[i] >= iMin && data[i] <= iMax){
//			msk[i] = 1;
//		}
//	}
//
//	// Guardar la mascara del indezx i en la plantilla de pixeles buenos tmp
//	tmp = tmp | (msk * (unsigned char)(1 << index));
//
//	// Emplear la mascara para cribar los pixeles malos
//	for(int i=0; i < size_data; i++){
//		data[i] = data[i] * (unsigned int)msk[i];
//	}
//}
//
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
//
//
//Mat getConstOlder(vector<Mat>& data, \
//             const Mat& tmp, \
//             Mat& pixCnt, \
//             int **disp) {
//
//	Mat con(data[0].size(), CV_64F);
//	vector<Mat> dat;
//
//	unsigned int loopCnt = 0;
//
//	// Calculo del logaritmo comun (base 10) de la imagen
//	dat.push_back(log10(data[0]));
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Calculo del logaritmo comun (base 10) de la imagen
//		dat.push_back(log10(data[iq]));
//
//		// Obtencion de la mascara
//		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			Mat mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//			int dx = disp[iq][0] - disp[ir][0];
//			int dy = disp[iq][1] - disp[ir][1];
//
//			// Calculo de los extremos de las ventanas
//			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
//			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
//			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
//			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS
//
//			// Calcular ventanas de mascara. MskiqROI y mskirROI son del mismo tama�o, aunque estan desplazadas unas con respecto a la otra una distancia relativa.
//			Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
//			Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Calcular la mascara de las ventanas
//			Mat mskDouble;
//			Mat msk = mskiqROI.mul(mskirROI);
//
//			msk.convertTo(mskDouble, CV_64F);
//
//			cout << "-------------------------------------" << endl;
//			cout << "Depth mskiq: " << mskiq.depth() << "\tDepth msk: " << msk.depth() << "\tDepth mskDouble: " << mskDouble.depth() << endl;
//
//			cout << "Depth dat[iq]: " << dat[iq].depth() << "\tDepth dat[ir]: " << dat[ir].depth() << endl;
//			//----------------------------------------------------------------------------------------
//			cout << "Pasada " << ir << endl << "-------------------------------------" << endl;
//
//			Size size = msk.size();
//			int xmax = size.height;
//			int ymax = size.width;
//
//			cout << "Msk \tAncho: " << xmax << "\tAlto: " << ymax << endl;
//			//----------------------------------------------------------------------------------------
//
//			// Calcular ventanas de datos
//			Mat datiqROI(dat[iq], Range(jyl, jyh), Range(jxl, jxh));
//			Mat datirROI(dat[ir], Range(iyl, iyh), Range(ixl, ixh));
//			//iq=log10(dataiq.at<uchar>(yActualQ, xActualQ))
//			//----------------------------------------------------------------------
//			size = datiqROI.size();
//			xmax = size.height;
//			ymax = size.width;
//			cout << "datiqROI\tAncho: " << xmax << "\tAlto: " << ymax << endl;
//
//			size = datirROI.size();
//			xmax = size.height;
//			ymax = size.width;
//			cout << "datirROI\tAncho: " << xmax << "\tAlto: " << ymax << endl;
//			//-----------------------------------------------------------------------
//
//
//			// Calcular diferencia
//			Mat diff = (datiqROI - datirROI).mul(mskDouble);
//
//			//-----------------------------------------------------------------------
//			size = diff.size();
//			xmax = size.height;
//			ymax = size.width;
//			cout << "diff\tAncho: " << xmax << "\tAlto: " << ymax << endl;
//			//-----------------------------------------------------------------------
///*
//			imshow("diff",diff);
//
//						waitKey( 0 );
//*/
////----------------------------------------------
//			double min,max;
//			Point min_loc,max_loc;
//
//			minMaxLoc(diff, &min, &max, &min_loc, &max_loc);
//
//			Mat B;
//
//			diff.convertTo(B, CV_32F);
//
//			diff.convertTo(B,CV_8U,255.0/(max-min),-255.0/min);
//
//			cout << "Min: " << min << "\tMax: " << max << endl;
//
////---------------------------------------------------
//
//			// Calcular ventanas del termino constante
//			Mat conJROI(con, Range(jyl, jyh), Range(jxl, jxh));
//			Mat conIROI(con, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Aplicar la diferencia a las ventanas del termino constante
//			conJROI = conJROI + diff;
//			conIROI = conIROI - diff;
//
//			// Calcular ventanas de la matriz de conteo de pares de pixeles
//			Mat pixCntJROI(pixCnt, Range(jyl, jyh), Range(jxl, jxh));
//			Mat pixCntIROI(pixCnt, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Aplicar la mascara a las ventanas de la matriz de pares de pixeles
//			pixCntJROI = pixCntJROI + msk;
//			pixCntIROI = pixCntIROI + msk;
//
//
//			//---------------------------------------------------------
//			double minVal, maxVal;
//			minMaxLoc(mskDouble, &minVal, &maxVal); //find minimum and maximum intensities
//
//			cout << "Depth CON: " << con.depth() << "\tDepth PixCnt: " << pixCnt.depth() << endl;
//			cout << "Channels CON: " << con.channels() << "\tChannels PixCnt: " << pixCnt.channels() << endl;
//			cout << "Min: " << minVal << "\tMax: " << maxVal << endl;
//			cout << "\nValor mskDouble: " << con.at<double>(Point(260,260)) << endl;
//
//			imshow("con",con);
//			waitKey( 0 );
//
//#ifdef PROGRESS
//
//			cout << "getConst: Iteracion " << ++loopCnt << " de 28..." << endl;
//
//#endif
//
//		}
//
//	}
//
//	///////////////////////////
//	int gh; cin >> gh;
//	///////////////////////////
//	return con;
//
//}
//
////*************************************************************************************
////*************************************************************************************
//void getCon(const Mat& mskiq,const Mat& mskir,const Mat& dataiq,const Mat& datair, Mat& pixCnt, Mat& con, int dx, int dy){
//
//	// Calculo de los extremos de las ventanas
//	unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
//	unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
//	unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
//	unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS
//
//	unsigned int maxY = (jyh - jyl), maxX = (jxh - jxl);		//Calculamos el tama�o de las ROI
//	unsigned int xActualQ, yActualQ, xActualR, yActualR;
//
//	for (unsigned int y = 0; y < maxY; y++){
//		for (unsigned int x = 0; x < maxX; x++){
//			//Calculamos los indices de desplazamientos actuales de ambas regiones
//			yActualQ = (y + jyl);	yActualR = (y + iyl);
//			xActualQ = (x + jxl);	xActualR = (x + ixl);
//
//			uchar msk = mskiq.at<uchar>(yActualQ, xActualQ) * mskir.at<uchar>(yActualR,xActualR);		//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!
//
//			//**************************Calculamos el CON**************************
//			double mskDouble = msk;
//
//			double logDataiQ = log10(dataiq.at<double>(yActualQ, xActualQ));
//			double logDataiR = log10(datair.at<double>(yActualR, xActualR));
//
//			//double diff = (logDataiQ - logDataiR) (AND o *) mskDouble;
//			double diff = (logDataiQ - logDataiR) * mskDouble;		//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!
//
//			con.at<double>(yActualQ, xActualQ) = con.at<double>(yActualQ, xActualQ) + diff;
//			con.at<double>(yActualR, xActualR) = con.at<double>(yActualR, xActualR) - diff;
//
//			//*************************Calculamos el PixCnt************************
//			pixCnt.at<uchar>(yActualQ, xActualQ) = pixCnt.at<uchar>(yActualQ, xActualQ) + msk;
//			pixCnt.at<uchar>(yActualR, xActualR) = pixCnt.at<uchar>(yActualR, xActualR) + msk;
//		}
//	}
//}
//
//
//
//Mat getConst(vector<Mat>& data, \
//             const Mat& tmp, \
//             Mat& pixCnt, \
//             int **disp) {
//
//	Mat con(data[0].size(), CV_64F);
//	vector<Mat> dat;
//
//	unsigned int loopCnt = 0;
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Obtencion de la mascara
//		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			Mat mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//			int dx = disp[iq][0] - disp[ir][0];
//			int dy = disp[iq][1] - disp[ir][1];
//
//			getCon(mskiq, mskir, data[iq], data[ir], pixCnt, con, dx, dy);
//
//			cout << "Pasada " << iq << "," << ir << endl;
//			imshow("con",con);
//
//			waitKey( 0 );
//
//#ifdef PROGRESS
//
//			cout << "getConst: Iteracion " << ++loopCnt << " de 28..." << endl;
//
//#endif
//
//		}
//
//	}
//
//	/////////////ELIMINAR/////////////
//	int gh; cin >> gh;
//	//////////////////////////////////
//	return con;
//
//}
//
////*************************************************************************************
////*************************************************************************************
//
//
//void doIterationNueva(const Mat& con, \
//                 Mat& gain, \
//                 const Mat& tmp, \
//                 const Mat& pixCnt, \
//                 const int disp[8][2]) {
//
//	unsigned int loopCnt = 0;
//
//	// Creacion de la ganancia temporal
//	Mat gainTmp;
//	con.copyTo(gainTmp);
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Obtencion de la mascara
//		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			Mat mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//			int dx = disp[iq][0] - disp[ir][0];
//			int dy = disp[iq][1] - disp[ir][1];
//
//			// Calculo de los extremos de las ventanas
//			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
//			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
//			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
//			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS
//
//			// Calcular ventanas de mascara
//			Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
//			Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Calcular la mascara de las ventanas
//			Mat msk = mskiqROI.mul(mskirROI);
//			msk.convertTo(msk, CV_64F);
//
//			// Calcular ventanas de ganancia y ganancia temporal
//			Mat gainTmpJROI(gainTmp, Range(jyl, jyh), Range(jxl, jxh));
//			Mat gainTmpIROI(gainTmp, Range(iyl, iyh), Range(ixl, ixh));
//			Mat gainJROI(gain, Range(jyl, jyh), Range(jxl, jxh));
//			Mat gainIROI(gain, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Modificar la ganancia temporal en base a la ganancia y la mascara
//			gainTmpJROI = gainTmpJROI + gainIROI.mul(msk);
//			gainTmpIROI = gainTmpIROI + gainJROI.mul(msk);
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
//	// Calcular ganancia unitaria para no dividir por cero
//	Mat pixCntAux = max(pixCnt, 1.0);
//	pixCntAux.convertTo(pixCntAux, CV_64F);
//
//	gainTmp = gainTmp / pixCntAux;
//
//	// Eliminar elementos a cero (de la matriz de pares de pixeles)
//	Mat index = min(pixCnt, 1.0);
//	index.convertTo(index, CV_64F);
//
//	gainTmp = gainTmp.mul(index);
//
//	// Calcular sumatorios
//	double sum2 = sum(gainTmp)[0];
//	double sum3 = sum(gainTmp.mul(gainTmp))[0];
//	double nPix = sum(index)[0];
//
//	// Eliminar elementos mas de 5-Sigma veces alejados de la media
//	double ave2 = sum2 / nPix; //Ganancia media de la CCD
//	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);
//
//	//Marco los elementos las posiciones de los elementos con un cero que son
//	index = (abs(gainTmp - ave2) > fiveSigma) / 255;
//	index.convertTo(index, CV_64F);
//
//	sum2 = sum2 - sum(gainTmp.mul(index))[0];
//	nPix = nPix - sum(index)[0];
//
//	// Normalizar la tabla de ganancias
//	ave2 = sum2 / nPix;
//
//	gainTmp = gainTmp - ave2;
//
//	// Devolver la tabla de ganancias
//	gain = gainTmp;
//
//}
//
//
////*************************************************************************************
////*************************************************************************************
//
//void getGainTmp(const Mat& mskiq,const Mat& mskir, Mat& gainTmp, Mat& gain, int dx, int dy){
//
//	// Calculo de los extremos de las ventanas
//	unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
//	unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
//	unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
//	unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS
//
//	unsigned int maxY = (jyh - jyl), maxX = (jxh - jxl);		//Calculamos el tama�o de las ROI
//	unsigned int xActualQ, yActualQ, xActualR, yActualR;
//
//	for (unsigned int y = 0; y < maxY; y++){
//		for (unsigned int x = 0; x < maxX; x++){
//			//Calculamos los indices de desplazamientos actuales de ambas regiones
//			yActualQ = (y + jyl);	yActualR = (y + iyl);
//			xActualQ = (x + jxl);	xActualR = (x + ixl);
//
//			//Obtenemos el valor de la mascara en el pixel actual
//
//			uchar msk = mskiq.at<uchar>(yActualQ, xActualQ) * mskir.at<uchar>(yActualR,xActualR);	//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!
//
//			double mskDouble = msk * REL8TO64;		//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!
//
//			//Obtenemos los valores de las ganancias en los pixeles actuales
//			double gainTmpJ = gainTmp.at<double>(yActualQ, xActualQ);
//			double gainTmpI = gainTmp.at<double>(yActualR, xActualR);
//			double gainJ = gain.at<double>(yActualQ, xActualQ);
//			double gainI = gain.at<double>(yActualR, xActualR);
//
//			//Modificar la ganancia temporal en base a la ganancia y la mascara
//			gainTmp.at<double>(yActualQ, xActualQ) = gainTmpJ + (gainI*mskDouble);
//			gainTmp.at<double>(yActualR, xActualR) = gainTmpI + (gainJ*mskDouble);
//
//			/*
//			//Modificar la ganancia temporal en base a la ganancia y la mascara DE OTRA FORMA
//			gainTmp.at<double>(yActualQ, xActualQ) = gainTmp.at<double>(yActualQ, xActualQ) + ( gain.at<double>(yActualR, xActualR) * mskDouble );
//			gainTmp.at<double>(yActualR, xActualR) = gainTmp.at<double>(yActualR, xActualR) + ( gain.at<double>(yActualQ, xActualQ) * mskDouble );
//			*/
//		}
//	}
//}
//
//void calculateStats(const Mat& pixCnt, Mat& gainTmp, Mat& gain){
//	Mat pixCntAux = max(pixCnt, 1.0);
//	pixCntAux.convertTo(pixCntAux, CV_64F);
//
//	Mat index = min(pixCnt, 1.0);
//	index.convertTo(index, CV_64F);
//
//	Size s = gainTmp.size();
//	int xmax = s.height;
//	int ymax = s.width;
//	for (int y = 0 ; y < ymax; y++)
//		for (int x = 0; x < xmax; x++)
//			gainTmp.at<double>(y, x) = gainTmp.at<double>(y, x) / pixCntAux.at<double>(y, x);
//
//	double sum2 = 0, sum3 = 0, nPix = 0;
//	for (int y = 0 ; y < ymax; y++){
//		for (int x = 0; x < xmax; x++){
//			sum2 += gainTmp.at<double>(y,x);												//COMPROBAR VAL[0]
//			sum3 += (gainTmp.at<double>(y,x) * gainTmp.at<double>(y,x));					//COMPROBAR VAL[0]
//			nPix += index.at<double>(y,x);													//COMPROBAR VAL[0]
//		}
//	}
//
//	// Eliminar elementos mas de 5-Sigma veces alejados de la media
//	double ave2 = sum2 / nPix; //Ganancia media de la CCD
//	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);
//
//	//Marco los elementos las posiciones de los elementos con un cero que son
//	for (int y = 0 ; y < ymax; y++){
//		for (int x = 0; x < xmax; x++){
//			index.at<double>(y,x) = (abs(gainTmp.at<double>(y,x) - ave2) > fiveSigma) / 255;
//		}
//	}
//	index.convertTo(index, CV_64F);
//
//	for (int y = 0 ; y < ymax; y++){
//		for (int x = 0; x < xmax; x++){
//			sum2 -= (gainTmp.at<double>(y,x) * index.at<double>(y,x));				//COMPROBAR VAL[0] EN AMBOS
//			nPix -= (index.at<double>(y,x));										//COMPROBAR VAL[0]
//		}
//	}
//
//	// Normalizar la tabla de ganancias
//	ave2 = sum2 / nPix;
//	for (int y = 0 ; y < ymax; y++){
//		for (int x = 0; x < xmax; x++){
//			gainTmp.at<double>(y,x) -= ave2;
//		}
//	}
//
//	// Devolver la tabla de ganancias
//	gain = gainTmp;
//}
//
//void doIteration(const Mat& con, \
//					  Mat& gain, \
//					  const Mat& tmp, \
//					  const Mat& pixCnt, \
//					  const int disp[8][2]) {
//
//	unsigned int loopCnt = 0;
//
//	// Creacion de la ganancia temporal
//	Mat gainTmp;
//	con.copyTo(gainTmp);
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Obtencion de la mascara
//		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			Mat mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//			int dx = disp[iq][0] - disp[ir][0];
//			int dy = disp[iq][1] - disp[ir][1];
//
//			getGainTmp(mskiq, mskir, gainTmp, gain, dx, dy);
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
//	calculateStats(pixCnt, gainTmp, gain);
//
//}
//
////*************************************************************************************
////*************************************************************************************
//
//
//Mat iterate(const Mat& con, \
//            Mat& gain, \
//            const Mat& tmp, \
//            const Mat& pixCnt, \
//            const int disp[8][2], \
//			const unsigned int loops) {
//
//	for(unsigned int i = 0; i < loops; i++) {
//
//		doIteration(con, gain, tmp, pixCnt, disp);
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
//	Mat flat = gain * log(10.0);
//	exp(flat, flat);
//
//	Mat tmpAux = (tmp > 0) / 255;
//	tmpAux.convertTo(tmpAux, CV_64F);
//
//	flat = flat.mul(tmpAux);
//
//	return flat;
//
//}



/*
int getImagesOriginal(vector<Mat>& data, \
              Mat& tmp, \
              const double rMin, \
              const double rMax) {

	// Comprobar el numero de imagenes
	unsigned int nImages = data.size();

	if(nImages > 8) {

		cerr << "Error - El numero de imagenes deben de ser ocho..." << endl;

		return 1;

	}

	// Calcular la plantilla de pixeles buenos
	for(unsigned int i = 0; i < 8; i++) {

		// Calcular la mascara, dividiendo por 255 debido a que, en OpenCV, TRUE
		// equivale a 255
		Mat msk = ((data[i] >= rMin) & (data[i] <= rMax)) / 255;

		cout << "Imagen " << i + 1 << ": " << (unsigned int)(sum(msk).val[0]) << " pixeles buenos..." << endl;

		// Guardar la mascara en la plantilla de pixeles buenos
		tmp = tmp | (msk * (1 << i));

		// Emplear la mascara para cribar los pixeles malos
		msk.convertTo(msk, data[i].type());
		data[i] = data[i].mul(msk);

	}

	return 0;

}
*/









































//
//char imageName[] = "./img/im0X.tiff";
//
//	vector<Mat> datacube;
//
//	// Leer imagenes desde fichero, guardandolas en el vector de datos
//	for(unsigned int i = 0; i < 8; i++) {
//
//		imageName[9] = 49 + i;
//
//		datacube.push_back(imread(imageName, -1));
//		datacube[i].convertTo(datacube[i], CV_64F);
//
//	}
//
//	// Tabla de desplazamientos - TODO Calcular los desplazamientos de forma automatica
//	const int disp[8][2] = {{0, 0}, {8, -26}, {32, -42}, {63, -30}, {79, -3}, {64, 24}, {35, 33}, {3, 18}};
//
//	// Declaracion de la plantilla de pixeles buenos
//	Mat tmp(datacube[0].size(), CV_8UC1, 0.0);
//
//	// Calculo de la plantilla de pixeles malos
//	int err = getImages(datacube, tmp, RMIN, RMAX);
//
//	if (err)
//		return 1;
//
//	// Declaracion de la matriz de conteo de pares de pixeles
//	Mat pixCnt(datacube[0].size(), CV_8UC1, 0.0);
//
//	// Calculo del termino constante
//	Mat con = getConst(datacube, tmp, pixCnt, disp);
//
//#ifdef DEBUG
//
//	namedWindow("Constant term", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
//	imshow("Constant term", con);
//	waitKey(0);
//
//#endif
//
//	// Calculo de la ganancia unitaria
//	Mat pixCntAux = max(pixCnt, 1.0);
//	pixCntAux.convertTo(pixCntAux, CV_64F);
//
//	Mat gain = con / pixCntAux;
//
//	// Calculo del flatfield
//	Mat flat = iterate(con, gain, tmp, pixCnt, disp, LOOPS);
//	flat = to16U(flat);
//
//#ifdef DEBUG
//
//	namedWindow("Flatfield", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
//	imshow("Flatfield", flat);
//	waitKey(0);
//
//#endif
//
//	imwrite("./out/ff.png", flat);
//
//	return 0;
//
//}
//
