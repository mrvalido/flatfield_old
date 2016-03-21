#include "utility.hpp"


ImageValDouble log_10(ImageValInt val){
	float ret;
	int size_val = val.size();

	ImageValDouble tmp(size_val);
	for(int i=0; i < size_val; i++){
		ret = log(val[i]);

		ret = max(ret,(float)0);

		ret /= log(10.0);
		tmp[i] = ret;
	}

	return tmp;
}

ImageValDouble Max(const ImageValDouble& val,double x){
	int size_val = val.size();
	ImageValDouble tmp(size_val);
	tmp=val;


	for(int i=0; i < size_val; i++){

		if (val[i] < x)
			tmp[i]=x;

	}
	return tmp;
}

ImageValDouble Min(const ImageValDouble& val,double x){
	int size_val = val.size();
	ImageValDouble tmp(size_val);
	tmp=val;


	for(int i=0; i < size_val; i++){

		if (val[i]  >= x)
			tmp[i]=x;
		else
			tmp[i]=val[i];

	}
	return tmp;
}


ImageValDouble toDouble(const ImageValChar& val) {

	int size_val = val.size();
	ImageValDouble ret(size_val);
	for(int i=0; i < size_val; i++){
		ret[i]=(double)val[i];
	}
	return ret;

}
//Mat to16U(const Mat& matrix) {
//
//	Mat ret;
//	matrix.copyTo(ret);
//
//	double vMin, vMax;
//
//	minMaxLoc(ret, &vMin, &vMax, NULL, NULL);
//	ret += abs(vMin);
//
//	minMaxLoc(ret, &vMin, &vMax, NULL, NULL);
//	ret /= vMax;
//	ret *= 65535;
//
//	ret.convertTo(ret, CV_16UC1);
//
//	return ret;
//
//}
