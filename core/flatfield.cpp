#include "flatfield.hpp"
#include "utility.hpp"


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
	unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + dimY; // FILAS
	unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + dimX; // COLUMNAS
//	unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + dimY; // FILAS
//	unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + dimX; // COLUMNAS

	cout << "jyl: " << jyl << "      jxl: " << jxl << "        jyh: " << jyh << "        jxh: " << jxh << endl;

	int ancho = (jxl-jxh);
	ImageValDouble ROI((jyl-jyh) * (jxl-jxh));

	// Calcular ventanas de mascara. MskiqROI y mskirROI son del mismo tamaño, aunque
	//estan desplazadas unas con respecto a la otra una distancia relativa.
	for(int y=jyh; y <= jyl; y++){
		for(int x=jxh; x <= jxl; x++){
			ROI[(y-jyh)*ancho + (x-jxh)] = (double)val[ind(y,x)];
		}
	}

	return ROI;
}

template <typename T>
void sumROI(valarray<T>& val, valarray<T>& ROI, int dx, int dy){

	// Calculo de los extremos de las ventanas
	unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + dimY; // FILAS
	unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + dimX; // COLUMNAS


	int ancho = (jxl-jxh);

	for(int y=jyh; y <= jyl; y++){
		for(int x=jxh; x <= jxl; x++){
			 val[ind(y,x)] = ROI[(y-jyh)*ancho + (x-jxh)];
		}
	}
}

ImageValDouble getConst(vector<ImageValInt>& data, const ImageValChar& tmp, ImageValDouble& pixCnt, int centros[8][2]) {

	vector<ImageValDouble> dat;
	ImageValDouble con(data[0].size());

	// Calculo del logaritmo comun (base 10) de la imagen
	dat.push_back(log_10(data[0]));

	for(unsigned int iq = 1; iq < 2; iq++) {

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
