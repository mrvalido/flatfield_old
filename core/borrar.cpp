ImageValDouble getConst(vector<ImageValInt>& data, const ImageValChar& tmp, ImageValDouble& pixCnt, int centros[8][2]) {

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



			Mat im(dimY, dimX, CV_8UC1, Scalar(0));  //Es un tipo de dato de 4 bytes 32S

				//Se pone primero el eje Y y despues el eje X
				for (long y=0; y<dimY; y++){
						for (long x=0; x<dimX; x++){
						im.at<uchar>(y,x) = (uchar)(mskir[ind( y, x )]*100);
					}
				}


//			ImageValChar im8 = escalado8(datiqROI);

//				Mat im(Alto, Ancho, CV_8UC1, Scalar(0));  //Es un tipo de dato de 4 bytes 32S
//
//				//Se pone primero el eje Y y despues el eje X
//				for (long y=0; y<Alto; y++){
//						for (long x=0; x<Ancho; x++){
//						im.at<uchar>(y,x) = (uchar)(mskir[y*Ancho+x]);
//					}
//				}


				imwrite("msk.jpeg", im);
				namedWindow( "Display window", WINDOW_NORMAL);// Create a window for display.
				imshow( "Display window", im );

				//waitKey(0);


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
