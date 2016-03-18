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

Mat to16U(const Mat& matrix) {

	Mat ret;
	matrix.copyTo(ret);

	double vMin, vMax;

	minMaxLoc(ret, &vMin, &vMax, NULL, NULL);
	ret += abs(vMin);

	minMaxLoc(ret, &vMin, &vMax, NULL, NULL);
	ret /= vMax;
	ret *= 65535;

	ret.convertTo(ret, CV_16UC1);

	return ret;

}
