#include"ImageProc.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>

Mat drawShadedSquare(TipoQuadrado tipo) {
	Mat square = Mat(256, 256, CV_8UC1);
	switch (tipo)
	{
	case TipoQuadrado::left: {  //square with shadow on left
		for (int i = 0; i < square.rows; i++) {
			for (int j = 0; j < square.cols; j++) {
				square.data[(j*square.cols) + i] = i;
			}
		}
		break;
	}
	case TipoQuadrado::top: {    //square with shadow on top
		for (int i = 0; i < square.rows; i++) {
			for (int j = 0; j < square.cols; j++) {
				square.data[(i*square.rows) + j] = i;
			}
		}
		break;
	}
	case corner: {         //square with shadow on corner
		for (int i = 0; i < square.rows; i++) {
			for (int j = 0; j < square.cols; j++) {
				square.data[(j*square.cols) + i] = ceil((i + j) / 2.0);
			}
		}
		break;
	}
	default:
		break;
	}

	return square;
}

Mat find_replace_value(Mat src,uint value,uint replace_with){

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		if (src.data[i] == value) {
			out.data[i] = replace_with;
		}else {
			out.data[i] = src.data[i];
		}
	}
	return out;
}

Mat vcpi_gray_negative(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[i] = 255 - src.data[i];
	}
	return out;
}

Mat vcpi_rgb_negative(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));

	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[i] = 255 - src.data[i];
		out.data[i + 1] = 255 - src.data[i + 1];
		out.data[i + 2] = 255 - src.data[i + 2];
	}

	return out;
}

Mat vcpi_rgb_remove_red(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));

	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[i] = src.data[i];
		out.data[i + 1] = src.data[i + 1];
		out.data[i + 2] = 0;
	}

	return out;
}

Mat vcpi_rgb_remove_blue(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));

	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[i] = 0;
		out.data[i + 1] = src.data[i + 1];
		out.data[i + 2] = src.data[i + 2];
	}

	return out;
}

Mat vcpi_rgb_remove_green(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));

	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[i] = src.data[i];
		out.data[i + 1] = 0;
		out.data[i + 2] = src.data[i + 2];
	}

	return out;
}

Mat vcpi_rgb_to_gray(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	uint j = 0;
	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[j] = src.data[i] * 0.114 + src.data[i + 1] * 0.587 + src.data[i + 2] * 0.299;
		j++;
	}

	return out;
}

Mat vcpi_rgb_get_R(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));
	uint j = 0;
	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[j] = src.data[i + 2];
		j++;
	}

	return out;
}

Mat vcpi_rgb_get_G(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));
	uint j = 0;
	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[j] = src.data[i + 1];
		j++;
	}

	return out;
}

Mat vcpi_rgb_get_B(Mat src) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));
	uint j = 0;
	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[j] = src.data[i];
		j++;
	}

	return out;
}

Mat vcpi_gray_to_binary(Mat src, int Lower_threshold, int Upper_threshold) {

	//check for input image
	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	//check for non-null variable inputs
	int Low_T = Lower_threshold == NULL ? 0 : Lower_threshold;
	int Up_T = Upper_threshold == NULL ? 0 : Upper_threshold;
	Mat out = Mat(src.rows, src.cols, CV_8UC1); //Output Binary Image

	for (uint i = 0; i < src.cols*src.rows; i ++) {
		if (src.data[i] >= Low_T && src.data[i] <= Up_T) {
			out.data[i] = 255;
		}else {
			out.data[i] = 0;
		}
	}

	return out;
}

Mat VCPI_Segmenta_Cor(Mat src, int Lower_h, int Upper_h, int Lower_s, int Upper_s,
	int Lower_v, int Upper_v) {

	//check for input image
	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	//check for non-null variable inputs
	int Low_H = Lower_h == NULL ? 0 : Lower_h;
	int Low_S = Lower_s == NULL ? 0 : Lower_s;
	int Low_V = Lower_v == NULL ? 0 : Lower_v;
	int Up_H = Upper_h == NULL ? 255 : Upper_h;
	int Up_S = Upper_s == NULL ? 255 : Upper_s;
	int Up_V = Upper_v == NULL ? 255 : Upper_v;

	//prompt the user of chosen HSV values in case of Null inputs
	cout << "Parameters H:" << Low_H << ", " << Up_H << ", S: " << Low_S << ", " << Up_S << ", V: " << Low_V << ", " << Up_V << endl;

	Mat im_HSV = Mat(src.rows, src.cols, CV_8UC3); //HSV image variable
	cvtColor(src, im_HSV, CV_BGR2HSV);             //convert RGB to HSV

	vector<Mat> HSV_Channels;                      //Vector for HSV channels
	split(im_HSV, HSV_Channels);                   //Split the HSV Image

	Mat im_HSV_thresholded = Mat(src.rows, src.cols, CV_8UC1);   //Thresholded HSV Image variable
	im_HSV_thresholded = vcpi_gray_to_binary(HSV_Channels[0],Low_H,Up_H) & vcpi_gray_to_binary(HSV_Channels[1], Low_S, Up_S) & vcpi_gray_to_binary(HSV_Channels[2], Low_V, Up_V); //Threshold
	
	Mat out = Mat(src.rows, src.cols, CV_8UC3);                  //RGB output Image

	//Mask the Threshold result with the input RGB IMage
	uint j = 0;
	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		out.data[i] = src.data[i] & im_HSV_thresholded.data[j];           //B
		out.data[i + 1] = src.data[i + 1] & im_HSV_thresholded.data[j];   //G
		out.data[i + 2] = src.data[i + 2] & im_HSV_thresholded.data[j];   //R
		j++;
	}
	return out;
}

Mat vcpi_scale_gray_to_rgb(Mat src){

	//check for input image
	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));
	uint j = 0;
	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()){
		//Blue COLOR
		if (src.data[i] >= 128){
			out.data[j] = 0;
		}else if (src.data[i] >= 64){
			out.data[j] = 255 - ((src.data[i] - 64) * 4); //rampa decrescente
		}else{
			out.data[j] = 255;
		}
		
		//Green COLOR
		if (src.data[i] >= 192){
			out.data[j + 1] = 255 - ((src.data[i] - 192) * 4); //rampa decrescente
		}else if (src.data[i] >= 64){
			out.data[j + 1] = 255;
		}else{
			out.data[j + 1] = src.data[i] * 4; //rampa crescente
		}
		
		//Red COLOR
		if (src.data[i] >= 192){
			out.data[j + 2] = 255;
		}else if (src.data[i] >= 128){
			out.data[j + 2] = (src.data[i] - 128) * 4; //rampa crescente
		}else{
			out.data[j + 2] = 0;
		}
		j+= out.channels();
	}
	return out;
}

Mat vcpi_convolucao(Mat src, Mat kernel) {


	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src; 
	}else if (kernel.empty()){          	//check for input kernel
		cout << "There is no kernel!" << endl;
		return src;
	}else if ((kernel.rows%2)==0 || (kernel.rows % 2) == 0)	{ 	//kernel must not have pair dimensions, but can have different dimensions
		cout << "Kernel dimensions must not be pair!" << endl;
		return src;
	}
	//find the kernel midpoint
	int kern_pad_x = (kernel.rows == 1 ? 0 : ((kernel.rows - 1) / 2));
	int kern_pad_y = (kernel.cols == 1 ? 0 : ((kernel.rows - 1) / 2));

	//get the sum of the kernel points
	double kernel_sum = 0.0f;
	for (uint i = 0; i < kernel.rows*kernel.cols; i++){
			kernel_sum += (double) kernel.data[i];
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

		for (uint x = kern_pad_x; x < (src.rows-kern_pad_x); x++) {           //run the image on one axis
			for (uint y = kern_pad_y; y < (src.cols- kern_pad_y); y++) {      //run the image on another axis

				double matrix_sum = 0.0f;                                              //variable to hold the matricial sum
				for (int kern_x = -kern_pad_x; kern_x <= kern_pad_x; kern_x++) {       //run the kernel on one axis
					for (int kern_y = -kern_pad_y; kern_y <= kern_pad_y; kern_y++){    //run the kernel on another axis
																					   //calculate the sum of points + kernel values surrounding the current point
						matrix_sum += src.at<uchar>(x + kern_x, y + kern_y)*kernel.at<uchar>(kern_x + kern_pad_x, kern_y + kern_pad_y); 
					}
				}
				out.at<uchar>(x, y) =(uchar) ((double) matrix_sum /(double) kernel_sum);
			}
		}

	return out;
}

Mat vcpi_median_filter(Mat src, uint kernel_size) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}	else if (kernel_size == 1) { 	//kernel must not be 1x1
		cout << "Kernel dimensions must be greater than 1 !" << endl;
		return src;
	}	else if ((kernel_size % 2) == 0) { 	//kernel must not have pair dimensions, but can have different dimensions
		cout << "Kernel dimensions must not be pair!" << endl;
		return src;
	}

	int kern_pad = ((kernel_size - 1) / 2);
	int array_size = kernel_size * kernel_size;

	uint8_t *temp = (uint8_t *)malloc(array_size * sizeof(uint8_t));

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
					temp[(kernel_size*(kern_x+kern_pad))+(kern_y+kern_pad)] = src.at<uchar>(x + kern_x, y + kern_y);		   //get the current point
				}
			}
			sort(temp,temp+array_size);
			out.at<uchar>(x, y) = (uchar) temp[kern_pad];
		}
	}

	free(temp);
	return out;
}

Mat vcpi_gray_to_binary_global_mean(Mat src) {

	//check for input image
	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}

	uint8_t threshold = 0;
	uint64_t acc = 0;
	for (uint i = 0; i < src.cols*src.rows; i++) {
		acc += src.data[i];
	}
	threshold = ((float)(acc)/(float)(src.cols*src.rows));

	Mat out = Mat(src.rows, src.cols, CV_8UC1); //Output Binary Image
	for (uint i = 0; i < src.cols*src.rows; i++) {
		if (src.data[i] >= threshold) {
			out.data[i] = 255;
		}
		else {
			out.data[i] = 0;
		}
	}
	return out;
}

Mat vcpi_gray_to_binary_midpoint(Mat src, uint kernel_size) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (kernel_size == 1) { 	//kernel must not be 1x1
		cout << "Kernel dimensions must be greater than 1 !" << endl;
		return src;
	}
	else if ((kernel_size % 2) == 0) { 	//kernel must not have pair dimensions, but can have different dimensions
		cout << "Kernel dimensions must not be pair!" << endl;
		return src;
	}

	int kern_pad = ((kernel_size - 1) / 2);
	int array_size = kernel_size * kernel_size;

	uint8_t *temp = (uint8_t *)malloc(array_size * sizeof(uint8_t));

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	uint8_t threshold = 0;
	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
					temp[(kernel_size*(kern_x + kern_pad)) + (kern_y + kern_pad)] = src.at<uchar>(x + kern_x, y + kern_y);		   //get the current point
				}
			}
			sort(temp, temp + array_size);
			threshold =  round((float)(temp[0]+temp[array_size-1])*0.5f);

			if (src.at<uchar>(x,y) >= threshold) {
				out.at<uchar>(x, y) = 255;
			}
			else {
				out.at<uchar>(x, y) = 0;
			}
		}
	}

	free(temp);
	return out;
}

Mat vcpi_gray_to_binary_bernsen(Mat src, uint kernel_size, uint Cmin) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (kernel_size == 1) { 	//kernel must not be 1x1
		cout << "Kernel dimensions must be greater than 1 !" << endl;
		return src;
	}
	else if ((kernel_size % 2) == 0) { 	//kernel must not have pair dimensions, but can have different dimensions
		cout << "Kernel dimensions must not be pair!" << endl;
		return src;
	}

	int kern_pad = ((kernel_size - 1) / 2);
	int array_size = kernel_size * kernel_size;

	uint8_t *temp = (uint8_t *)malloc(array_size * sizeof(uint8_t));

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	uint8_t threshold = 0;
	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
					temp[(kernel_size*(kern_x + kern_pad)) + (kern_y + kern_pad)] = src.at<uchar>(x + kern_x, y + kern_y);		   //get the current point
				}
			}
			sort(temp, temp + array_size);
			if ((temp[array_size - 1] - temp[0]) < Cmin) {
				threshold = 128;
			}
			else {
				threshold = round((float)(temp[0] + temp[array_size - 1])*0.5f);
			}
			

			if (src.at<uchar>(x, y) >= threshold) {
				out.at<uchar>(x, y) = 255;
			}
			else {
				out.at<uchar>(x, y) = 0;
			}
		}
	}

	free(temp);
	return out;
}

Mat vcpi_gray_to_binary_niblack(Mat src, uint kernel_size, float k) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (kernel_size == 1) { 	//kernel must not be 1x1
		cout << "Kernel dimensions must be greater than 1 !" << endl;
		return src;
	}
	else if ((kernel_size % 2) == 0) { 	//kernel must not have pair dimensions, but can have different dimensions
		cout << "Kernel dimensions must not be pair!" << endl;
		return src;
	}

	int kern_pad = ((kernel_size - 1) / 2);
	int array_size = kernel_size * kernel_size;

	uint8_t *temp = (uint8_t *)malloc(array_size * sizeof(uint8_t));

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	uint8_t threshold = 0;
	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
					temp[(kernel_size*(kern_x + kern_pad)) + (kern_y + kern_pad)] = src.at<uchar>(x + kern_x, y + kern_y);		   //get the current point
				}
			}
			uint64_t acc_mean = 0;
			for (uint i = 0; i < array_size; i++) {
				acc_mean += temp[i];
			}
			float mean = ((float)acc_mean) / ((float)(array_size));
			
			long double acc_std = 0.0;
			for (uint i = 0; i < array_size; i++) {
				acc_std = ((double)(temp[i]-mean))*((double)(temp[i] - mean));
			}
			double std = sqrt(((double)acc_std) / ((double)array_size));

			threshold = (uint8_t) round((mean + ((double)std*k)));

			if (src.at<uchar>(x, y) >= threshold) {
				out.at<uchar>(x, y) = 255;
			}
			else {
				out.at<uchar>(x, y) = 0;
			}
		}
	}

	free(temp);
	return out;
}