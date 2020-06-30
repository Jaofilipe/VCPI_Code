#include"ImageProc.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <deque>
#include <string> 

using namespace cv;
using namespace std;

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
	case TipoQuadrado::corner: {         //square with shadow on corner
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

Mat find_replace_value(Mat src, uint value, uint replace_with) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		if (src.data[i] == value) {
			out.data[i] = replace_with;
		}
		else {
			out.data[i] = src.data[i];
		}
	}
	return out;
}

Mat find_create_with(Mat src, uint value, uint red, uint green, uint blue) {

	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));
	uint helper = 0;
	for (uint i = 0; i < src.cols * src.rows ; i++) {
		if (src.data[i] == value) {
			out.data[helper] = blue ;
			out.data[helper+1] = green;
			out.data[helper+2] = red;
		}
		else {
			out.data[helper] = 0;
			out.data[helper + 1] = 0;
			out.data[helper + 2] = 0;
		}
		helper += out.channels();
	}
	return out;
}

Mat find_apply_created(Mat src, Mat src_created) {

	if (src.empty() || src_created.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));
	
	for (uint i = 0; i < src.cols * src.rows * src.channels(); i++) {
		out.data[i] = (src_created.data[i] > 0) ? src_created.data[i] : src.data[i];
			
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

	for (uint i = 0; i < src.cols*src.rows; i++) {
		if (src.data[i] >= Low_T && src.data[i] <= Up_T) {
			out.data[i] = 255;
		}
		else {
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
	//cout << "Parameters H:" << Low_H << ", " << Up_H << ", S: " << Low_S << ", " << Up_S << ", V: " << Low_V << ", " << Up_V << endl;

	Mat im_HSV = Mat(src.rows, src.cols, CV_8UC3); //HSV image variable
	cvtColor(src, im_HSV, CV_BGR2HSV);             //convert RGB to HSV

	vector<Mat> HSV_Channels;                      //Vector for HSV channels
	split(im_HSV, HSV_Channels);                   //Split the HSV Image

	Mat im_HSV_thresholded = Mat(src.rows, src.cols, CV_8UC1);   //Thresholded HSV Image variable
	im_HSV_thresholded = vcpi_gray_to_binary(HSV_Channels[0], Low_H, Up_H) & vcpi_gray_to_binary(HSV_Channels[1], Low_S, Up_S) & vcpi_gray_to_binary(HSV_Channels[2], Low_V, Up_V); //Threshold

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

Mat vcpi_scale_gray_to_rgb(Mat src) {

	//check for input image
	if (src.empty()) {
		cout << "There is no image!" << endl;
		return src;
	}
	Mat out = Mat(src.rows, src.cols, CV_8UC3, Scalar(0));
	uint j = 0;
	for (uint i = 0; i < src.cols*src.rows*src.channels(); i += src.channels()) {
		//Blue COLOR
		if (src.data[i] >= 128) {
			out.data[j] = 0;
		}
		else if (src.data[i] >= 64) {
			out.data[j] = 255 - ((src.data[i] - 64) * 4); //rampa decrescente
		}
		else {
			out.data[j] = 255;
		}

		//Green COLOR
		if (src.data[i] >= 192) {
			out.data[j + 1] = 255 - ((src.data[i] - 192) * 4); //rampa decrescente
		}
		else if (src.data[i] >= 64) {
			out.data[j + 1] = 255;
		}
		else {
			out.data[j + 1] = src.data[i] * 4; //rampa crescente
		}

		//Red COLOR
		if (src.data[i] >= 192) {
			out.data[j + 2] = 255;
		}
		else if (src.data[i] >= 128) {
			out.data[j + 2] = (src.data[i] - 128) * 4; //rampa crescente
		}
		else {
			out.data[j + 2] = 0;
		}
		j += out.channels();
	}
	return out;
}

Mat vcpi_convolucao(Mat src, Mat kernel) {


	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (kernel.empty()) {          	//check for input kernel
		cout << "There is no kernel!" << endl;
		return src;
	}
	else if ((kernel.rows % 2) == 0 || (kernel.rows % 2) == 0) { 	//kernel must not have pair dimensions, but can have different dimensions
		cout << "Kernel dimensions must not be pair!" << endl;
		return src;
	}
	//find the kernel midpoint
	int kern_pad_x = (kernel.rows == 1 ? 0 : ((kernel.rows - 1) / 2));
	int kern_pad_y = (kernel.cols == 1 ? 0 : ((kernel.rows - 1) / 2));

	//get the sum of the kernel points
	double kernel_sum = 0.0f;
	for (uint i = 0; i < kernel.rows*kernel.cols; i++) {
		kernel_sum += (double)kernel.data[i];
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	for (uint x = kern_pad_x; x < (src.rows - kern_pad_x); x++) {           //run the image on one axis
		for (uint y = kern_pad_y; y < (src.cols - kern_pad_y); y++) {      //run the image on another axis

			double matrix_sum = 0.0f;                                              //variable to hold the matricial sum
			for (int kern_x = -kern_pad_x; kern_x <= kern_pad_x; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad_y; kern_y <= kern_pad_y; kern_y++) {    //run the kernel on another axis
																				   //calculate the sum of points + kernel values surrounding the current point
					matrix_sum += src.at<uchar>(x + kern_x, y + kern_y)*kernel.at<uchar>(kern_x + kern_pad_x, kern_y + kern_pad_y);
				}
			}
			out.at<uchar>(x, y) = (uchar)((double)matrix_sum / (double)kernel_sum);
		}
	}

	return out;
}

Mat vcpi_convolucao_freqs(Mat src, Mat kernel, double divide_by) {
	 
	//same as convolucao but supports negative elements on matrix

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (kernel.empty()) {          	//check for input kernel
		cout << "There is no kernel!" << endl;
		return src;
	}
	else if ((kernel.rows % 2) == 0 || (kernel.rows % 2) == 0) { 	//kernel must not have pair dimensions, but can have different dimensions
		cout << "Kernel dimensions must not be pair!" << endl;
		return src;
	}
	//find the kernel midpoint
	int kern_pad_x = (kernel.rows == 1 ? 0 : ((kernel.rows - 1) / 2));
	int kern_pad_y = (kernel.cols == 1 ? 0 : ((kernel.rows - 1) / 2));

	Mat temp = Mat(src.rows, src.cols, CV_32FC1, Scalar(0));

	for (uint x = kern_pad_x; x < (src.rows - kern_pad_x); x++) {           //run the image on one axis
		for (uint y = kern_pad_y; y < (src.cols - kern_pad_y); y++) {      //run the image on another axis

			double matrix_sum = 0.0f;                                              //variable to hold the matricial sum
			for (int kern_x = -kern_pad_x; kern_x <= kern_pad_x; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad_y; kern_y <= kern_pad_y; kern_y++) {    //run the kernel on another axis
																				   //calculate the sum of points + kernel values surrounding the current point
					matrix_sum += src.at<uchar>(x + kern_x, y + kern_y) * kernel.at<schar>(kern_x + kern_pad_x, kern_y + kern_pad_y);
				}
			}
			temp.at<float>(x, y) = (float)((double)matrix_sum / (double)divide_by);
		}
	}
	uint max_value = 0;  //find the maximum value of function to equalize the data 
	uint holder = 0;
	for (uint i = 0; i < src.cols * src.rows; i++) {
		holder = (uint)round(temp.data[i]);
		max_value = (holder > max_value) ? holder : max_value;
	}


	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	if (max_value > 255) {
		for (uint i = 0; i < src.cols * src.rows; i++) {
			out.data[i] = (temp.data[i] > 0) ? ((uint)round(((temp.data[i] * 255.0) / ((double)max_value)))) : ((uint)0);
		}
	}
	else {
		for (uint i = 0; i < src.cols * src.rows; i++) {
			out.data[i] = (temp.data[i] > 0) ? ((uint)round(temp.data[i])) : ((uint)0);
		}

	}


	
return out;
}

Mat vcpi_median_filter(Mat src, uint kernel_size) {

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

	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
					temp[(kernel_size*(kern_x + kern_pad)) + (kern_y + kern_pad)] = src.at<uchar>(x + kern_x, y + kern_y);		   //get the current point
				}
			}
			sort(temp, temp + array_size);
			out.at<uchar>(x, y) = (uchar)temp[kern_pad];
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
	threshold = ((float)(acc) / (float)(src.cols*src.rows));

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
			threshold = round((float)(temp[0] + temp[array_size - 1])*0.5f);

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
				acc_std = ((double)(temp[i] - mean))*((double)(temp[i] - mean));
			}
			double std = sqrt(((double)acc_std) / ((double)array_size));

			threshold = (uint8_t)round((mean + ((double)std*k)));

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

Mat vcpi_gray_to_binary_Region_Growing(Mat src, uint x, uint y, uint Lower_Limit, uint Upper_Limit, TipoVizinhanca neighborhood) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if ((neighborhood != Quatro) && (neighborhood != Oito)) {
		cout << "Please select a valid pixel neighborhood" << endl;
		return src;
	}
	else if ((x > src.rows) || (y >= src.cols)) {
		cout << "Please select a valid pixel position" << endl;
		return src;
	}

	uint8_t unchecked = 127;

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(unchecked));

	uint8_t seed = src.at<uchar>(x, y);

	std::deque<coordinates> fifo;
	coordinates next_point;
	next_point.x = x;
	next_point.y = y;
	fifo.push_front(next_point);

	coordinates cur_point;

	while (!fifo.empty()) {

		cur_point = fifo.back();
		fifo.pop_back();

		if ((src.at<uchar>(cur_point.x, cur_point.y) <= (seed + Upper_Limit)) && (src.at<uchar>(cur_point.x, cur_point.y) >= (seed - Lower_Limit)) && (out.at<uchar>(cur_point.x, cur_point.y) == unchecked)) {

			out.at<uchar>(cur_point.x, cur_point.y) = 255;
			//cout << "x:"<< (int)cur_point.x << endl;
			//cout << "y:"<< (int)cur_point.y << endl;
			//cout << "\n" << endl;
			if (cur_point.x == 0) {
				if (cur_point.y == 0) {

					next_point.x = cur_point.x + 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}
					next_point.x = cur_point.x;
					next_point.y = cur_point.y + 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}
				}
				else if (cur_point.y == (src.cols - 1)) {

					next_point.x = cur_point.x + 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y - 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}
				else {

					next_point.x = cur_point.x + 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y - 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y + 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}

			}
			else if (cur_point.x == (src.rows - 1)) {
				if (cur_point.y == 0) {

					next_point.x = cur_point.x - 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y + 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}
				else if (cur_point.y == (src.cols - 1)) {

					next_point.x = cur_point.x - 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y - 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}
				else {

					next_point.x = cur_point.x - 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y + 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y - 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}
			}
			else {
				if (cur_point.y == 0) {

					next_point.x = cur_point.x + 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x - 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y + 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}
				else if (cur_point.y == (src.cols - 1)) {

					next_point.x = cur_point.x + 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x - 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y - 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}
				else {

					next_point.x = cur_point.x + 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x - 1;
					next_point.y = cur_point.y;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y + 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

					next_point.x = cur_point.x;
					next_point.y = cur_point.y - 1;
					if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
						fifo.push_front(next_point);
					}

				}
			}

			if (neighborhood == Oito) {
				if (cur_point.x == 0) {
					if (cur_point.y == 0) {

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}
					}
					else if (cur_point.y == (src.cols - 1)) {

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

					}
					else {

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}
					}
				}
				else if (cur_point.x == (src.rows - 1)) {
					if (cur_point.y == 0) {

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}
					}
					else if (cur_point.y == (src.cols - 1)) {

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}
					}
					else {

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}
					}
				}
				else {
					if (cur_point.y == 0) {

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}


					}
					else if (cur_point.y == (src.cols - 1)) {

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}
					}
					else {

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y + 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

						next_point.x = cur_point.x + 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

						next_point.x = cur_point.x - 1;
						next_point.y = cur_point.y - 1;
						if (out.at<uchar>(next_point.x, next_point.y) == unchecked) {
							fifo.push_front(next_point);
						}

					}
				}
			}
		}
		else {
			if (out.at<uchar>(cur_point.x, cur_point.y) == unchecked) { // if the pixel was unchecked
				out.at<uchar>(cur_point.x, cur_point.y) = 0;
			}
		}
	}
	for (uint i = 0; i < out.cols*out.rows; i++) {
		if (out.data[i] == unchecked) {
			out.data[i] = 0;
		}
	}
	return out;
}

Mat vcpi_gray_to_binary_otsu(Mat src) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}

	uint64_t gray_values[256] = { 0 };
	long double Within_class_Variance[256] = { 0.0f };

	uint64_t img_size = src.cols * src.rows;

	//get the histogram of the image
	for (uint i = 0; i < img_size; i++) {
		gray_values[src.data[i]]++;
	}

	//run through the histogram of the image
	for (uint i = 0; i <= 255; i++) {

		//for background
		//calcular o weight
		long double weight_background = 0.0f;
		for (uint u = 0; u < i; u++) {
			weight_background += gray_values[u];
		}
		weight_background = (double)((double)weight_background) / ((double)img_size);

		//calculate mean
		double mean_background = 0.0f;
		uint64_t holder = 0;  //set holder variable
		for (uint u = 0; u < i; u++) {
			mean_background += (u * gray_values[u]);
			holder += gray_values[u];
		}
		mean_background = (holder == 0 ? 0 : (double)((double)mean_background) / ((double)holder));

		//calculate variance
		long double variance_background = 0.0;
		for (uint u = 0; u < i; u++) {
			variance_background += pow(((double)(u - mean_background)), 2) * ((double)gray_values[u]);
		}

		variance_background = (holder == 0 ? 0 : ((double)variance_background) / ((double)holder));

		//for foreground
		//calcular o weight
		long double weight_foreground = 0.0f;
		for (uint u = i; u <= 255; u++) {
			weight_foreground += gray_values[u];
		}
		weight_foreground = (double)((double)weight_foreground) / ((double)img_size);
		//calculate mean
		double mean_foreground = 0.0f;
		holder = 0; //reset holder variable
		for (uint u = i; u <= 255; u++) {
			mean_foreground += (u * gray_values[u]);
			holder += gray_values[u];
		}
		mean_foreground = (holder == 0 ? 0 : (double)((double)mean_foreground) / ((double)holder));
		//calculate variance
		long double variance_foreground = 0.0;
		for (uint u = i; u <= 255; u++) {
			variance_foreground += pow(((double)(u - mean_foreground)), 2) * ((double)gray_values[u]);
		}
		variance_foreground = (holder == 0 ? 0 : ((double)variance_foreground) / ((double)holder));

		Within_class_Variance[i] = (long double)((weight_background * variance_background) + (weight_foreground * variance_foreground));
	}

	double min = Within_class_Variance[0];
	uint threshold = 0;
	for (uint u = 0; u <= 255; u++) {

		if (min > Within_class_Variance[u]) {
			min = Within_class_Variance[u];
			threshold = u;
		}
	}
	//cout << "Threshold:   " << threshold << endl;
	return vcpi_gray_to_binary(src, threshold);
}

Mat vcpi_dilation(Mat src, TipoVizinhanca neighborhood) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if ((neighborhood != Quatro) && (neighborhood != Oito)) {
		cout << "Please select a valid pixel neighborhood" << endl;
		return src;
	}
	uint8_t array_size = (uint8_t)neighborhood + 1;

	uint8_t *temp = (uint8_t *)malloc(array_size * sizeof(uint8_t));

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	if (neighborhood == Quatro) {

		for (uint x = 0; x < src.rows; x++) {
			for (uint y = 0; y < src.cols; y++) {

				fill(temp, temp + array_size, 0); //zero the array
				if (x == 0) {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);
					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else {
					
						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}

				}
				else if (x == (src.rows - 1)) {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
				}
				else {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						temp[4] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
				}
			}

		}
	}
	else if (neighborhood == Oito) {

		for (uint x = 0; x < src.rows; x++) {
			for (uint y = 0; y < src.cols; y++) {

				fill(temp, temp + array_size, 0); //set the array

				if (x == 0) {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						temp[3] = src.at<uchar>(x + 1, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x + 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						temp[4] = src.at<uchar>(x + 1, y + 1);
						temp[5] = src.at<uchar>(x + 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
				}
				else if (x == (src.rows - 1)) {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						temp[3] = src.at<uchar>(x - 1, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						temp[4] = src.at<uchar>(x - 1, y + 1);
						temp[5] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
				}
				else {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y + 1);
						temp[4] = src.at<uchar>(x + 1, y + 1);
						temp[5] = src.at<uchar>(x - 1, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						temp[4] = src.at<uchar>(x + 1, y - 1);
						temp[5] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						temp[4] = src.at<uchar>(x, y + 1);
						temp[5] = src.at<uchar>(x + 1, y + 1);
						temp[6] = src.at<uchar>(x - 1, y + 1);
						temp[7] = src.at<uchar>(x + 1, y - 1);
						temp[8] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*max_element(temp, temp + array_size);

					}
				}
			}
		}
	}


free(temp);
return out;
}

Mat vcpi_erosion(Mat src, TipoVizinhanca neighborhood){
	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if ((neighborhood != Quatro) && (neighborhood != Oito)) {
		cout << "Please select a valid pixel neighborhood" << endl;
		return src;
	}
	uint8_t array_size = (uint)neighborhood + 1;

	uint8_t *temp = (uint8_t *)malloc(array_size * sizeof(uint8_t));
	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	if (neighborhood == Quatro) {

		for (uint x = 0; x < src.rows; x++) {
			for (uint y = 0; y < src.cols; y++) {

				fill(temp, temp + array_size, 255); //set the array

				if (x == 0) {
					if (y == 0) {
						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);
					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}

				}
				else if (x == (src.rows - 1)) {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
				}
				else {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						temp[4] = src.at<uchar>(x, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
				}
			}

		}
	}
	else if (neighborhood == Oito) {

		for (uint x = 0; x < src.rows; x++) {
			for (uint y = 0; y < src.cols; y++) {

				fill(temp, temp + array_size, 255); //set the array

				if (x == 0) {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						temp[3] = src.at<uchar>(x + 1, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x + 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						temp[4] = src.at<uchar>(x + 1, y + 1);
						temp[5] = src.at<uchar>(x + 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
				}
				else if (x == (src.rows - 1)) {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y + 1);
						temp[3] = src.at<uchar>(x - 1, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x - 1, y);
						temp[2] = src.at<uchar>(x, y - 1);
						temp[3] = src.at<uchar>(x, y + 1);
						temp[4] = src.at<uchar>(x - 1, y + 1);
						temp[5] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
				}
				else {
					if (y == 0) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y + 1);
						temp[4] = src.at<uchar>(x + 1, y + 1);
						temp[5] = src.at<uchar>(x - 1, y + 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else if (y == (src.cols - 1)) {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						temp[4] = src.at<uchar>(x + 1, y - 1);
						temp[5] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
					else {

						temp[0] = src.at<uchar>(x, y);
						temp[1] = src.at<uchar>(x + 1, y);
						temp[2] = src.at<uchar>(x - 1, y);
						temp[3] = src.at<uchar>(x, y - 1);
						temp[4] = src.at<uchar>(x, y + 1);
						temp[5] = src.at<uchar>(x + 1, y + 1);
						temp[6] = src.at<uchar>(x - 1, y + 1);
						temp[7] = src.at<uchar>(x + 1, y - 1);
						temp[8] = src.at<uchar>(x - 1, y - 1);
						out.at<uchar>(x, y) = (uint8_t)*min_element(temp, temp + array_size);

					}
				}
			}
		}
	}

	free(temp);
	return out;
}

Mat vcpi_binary_close(Mat src, uint dilation_times, TipoVizinhanca dilation_neighborhood, uint erosion_times, TipoVizinhanca erosion_neighborhood) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if ((dilation_neighborhood != Quatro) && (dilation_neighborhood != Oito)) {
		cout << "Please select a valid pixel dilation neighborhood" << endl;
		return src;
	}
	else if ((erosion_neighborhood != Quatro) && (erosion_neighborhood != Oito)) {
		cout << "Please select a valid pixel erosion neighborhood" << endl;
		return src;
	}
	else if ((dilation_times<0) || (erosion_times<0)) {
		cout << "Please select a valid pixel erosion and dilation value" << endl;
		return src;
	}
	if (dilation_times == 0) {
		cout << "Warning, number of dilations is zero" << endl;
	}
	else if (erosion_times == 0) {
		cout << "Warning, number of erosions is zero" << endl;
	}

	Mat temp_one = src.clone();//since 0 is supported as erosion and dilation this avoids outputting the wrong matrix
	Mat temp_two = src.clone();

	for (uint i = 0; i < dilation_times; i++){     // this cycle improves execution time by avoiding to zero an array and performing a simple calculation
		if (i%2 == 0){
			temp_one = vcpi_dilation(temp_two, dilation_neighborhood);

		}else{
			temp_two = vcpi_dilation(temp_one, dilation_neighborhood);

		}
		
	}
	if (dilation_times %2 == 1){ //since 0 is supported as erosion and dilation this avoids outputting the wrong matrix
		temp_two = temp_one;
	}
	else {
		temp_one = temp_two;
	}

	for (uint i = 0; i < erosion_times; i++)
	{
		if (i % 2 == 0) {
			temp_one = vcpi_erosion(temp_two, erosion_neighborhood);

		}
		else {
			temp_two = vcpi_erosion(temp_one, erosion_neighborhood);

		}
	}

	return ((erosion_times % 2) == 1 ? temp_one : temp_two); //choose the correct image to output
}

Mat vcpi_binary_open(Mat src, uint erosion_times, TipoVizinhanca erosion_neighborhood, uint dilation_times, TipoVizinhanca dilation_neighborhood) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if ((dilation_neighborhood != Quatro) && (dilation_neighborhood != Oito)) {
		cout << "Please select a valid pixel dilation neighborhood" << endl;
		return src;
	}
	else if ((erosion_neighborhood != Quatro) && (erosion_neighborhood != Oito)) {
		cout << "Please select a valid pixel erosion neighborhood" << endl;
		return src;
	}
	else if ((dilation_times < 0) || (erosion_times < 0)) {
		cout << "Please select a valid pixel erosion and dilation value" << endl;
		return src;
	}
	if (dilation_times == 0) {
		cout << "Warning, number of dilations is zero" << endl;
	}
	else if (erosion_times == 0) {
		cout << "Warning, number of erosions is zero" << endl;
	}

	Mat temp_one = src.clone();  //since 0 is supported as erosion and dilation this avoids outputting the wrong matrix
	Mat temp_two = src.clone();

	for (uint i = 0; i < erosion_times; i++) {
		if (i % 2 == 0) {
			temp_one = vcpi_erosion(temp_two, erosion_neighborhood);

		}
		else {
			temp_two = vcpi_erosion(temp_one, erosion_neighborhood);

		}
	}
	
	if (erosion_times % 2 == 1){  //since 0 is supported as erosion and dilation this avoids outputting the wrong matrix
		temp_two = temp_one;
	}
	else {
		temp_one = temp_two;
	}

	for (uint i = 0; i < dilation_times; i++)
	{
		if (i % 2 == 0) {     // this cycle improves execution time by avoiding to zero an array and performing a simple calculation
			temp_one = vcpi_dilation(temp_two, dilation_neighborhood);

		}
		else {
			temp_two = vcpi_dilation(temp_one, dilation_neighborhood);

		}
	}

	return ((dilation_times % 2) == 1 ? temp_one : temp_two); //choose the correct image to output
}

Mat vcpi_gray_edge_prewitt(Mat src, float th) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	
	float percent_threshold = th == NULL ? 0.5f : th;

	Mat prewitt_vertical(3, 3, CV_8SC1, Scalar(0));
	prewitt_vertical.at<schar>(0, 0) = -1;
	prewitt_vertical.at<schar>(0, 1) = -1;
	prewitt_vertical.at<schar>(0, 2) = -1;
	prewitt_vertical.at<schar>(2, 0) = 1;
	prewitt_vertical.at<schar>(2, 1) = 1;
	prewitt_vertical.at<schar>(2, 2) = 1;

	Mat prewitt_horizontal(3, 3, CV_8SC1, Scalar(0));
	prewitt_horizontal.at<schar>(0, 0) = -1;
	prewitt_horizontal.at<schar>(1, 0) = -1;
	prewitt_horizontal.at<schar>(2, 0) = -1;
	prewitt_horizontal.at<schar>(0, 2) = 1;
	prewitt_horizontal.at<schar>(1, 2) = 1;
	prewitt_horizontal.at<schar>(2, 2) = 1;

	uint kernel_size = prewitt_horizontal.cols;
	int kern_pad = ((kernel_size - 1) / 2);
	int array_size = kernel_size * kernel_size;

	Mat vertical = Mat(src.rows, src.cols, CV_32S, Scalar(0));
	Mat horizontal = Mat(src.rows, src.cols, CV_32S, Scalar(0));

	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			int32_t matrix_sum = 0;  //variable to hold the matricial sum

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
																				   //calculate the sum of points + kernel values surrounding the current point
					matrix_sum += src.at<uchar>(x + kern_x, y + kern_y)*prewitt_vertical.at<schar>(kern_x + kern_pad, kern_y + kern_pad);
				}
			}
			vertical.at<int>(x, y) = (int)matrix_sum;
		}
	}

	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			int32_t matrix_sum = 0;  //variable to hold the matricial sum

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
																				   //calculate the sum of points + kernel values surrounding the current point
					matrix_sum += src.at<uchar>(x + kern_x, y + kern_y)*prewitt_horizontal.at<schar>(kern_x + kern_pad, kern_y + kern_pad);
				}
			}
			horizontal.at<int>(x, y) = (int)matrix_sum;
		}
	}

	Mat temp = Mat(src.rows, src.cols, CV_32S, Scalar(0));
	uint32_t max = 0;
	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			temp.at<int>(x, y) = (int32_t)sqrt(pow(vertical.at<int>(x, y), 2) + (pow(horizontal.at<int>(x, y), 2)));
			max = (max >= temp.at<int>(x, y) ? max : temp.at<int>(x, y));

		}
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));
	uint threshold = round(max*percent_threshold);

	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			if (temp.at<int>(x, y) >= threshold) {
				out.at<uchar>(x, y) = 255;
			}
			else {
				out.at<uchar>(x, y) = 0;
			}
		}
	}

	return out;
}

Mat vcpi_gray_edge_sobel(Mat src, float th) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}

	float percent_threshold = th == NULL ? 0.5f : th;

	Mat sobel_vertical(3, 3, CV_8SC1, Scalar(0));
	sobel_vertical.at<schar>(0, 0) = -1;
	sobel_vertical.at<schar>(0, 1) = -2;
	sobel_vertical.at<schar>(0, 2) = -1;
	sobel_vertical.at<schar>(2, 0) = 1;
	sobel_vertical.at<schar>(2, 1) = 2;
	sobel_vertical.at<schar>(2, 2) = 1;

	Mat sobel_horizontal(3, 3, CV_8SC1, Scalar(0));
	sobel_horizontal.at<schar>(0, 0) = -1;
	sobel_horizontal.at<schar>(1, 0) = -2;
	sobel_horizontal.at<schar>(2, 0) = -1;
	sobel_horizontal.at<schar>(0, 2) = 1;
	sobel_horizontal.at<schar>(1, 2) = 2;
	sobel_horizontal.at<schar>(2, 2) = 1;

	uint kernel_size = sobel_horizontal.cols;
	int kern_pad = ((kernel_size - 1) / 2);
	int array_size = kernel_size * kernel_size;

	Mat vertical = Mat(src.rows, src.cols, CV_32S, Scalar(0));
	Mat horizontal = Mat(src.rows, src.cols, CV_32S, Scalar(0));

	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			int32_t matrix_sum = 0;  //variable to hold the matricial sum

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
																				   //calculate the sum of points + kernel values surrounding the current point
					matrix_sum += src.at<uchar>(x + kern_x, y + kern_y)*sobel_vertical.at<schar>(kern_x + kern_pad, kern_y + kern_pad);
				}
			}
			vertical.at<int>(x,y) = (int)matrix_sum;
		}
	}

	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			int32_t matrix_sum = 0;  //variable to hold the matricial sum

			for (int kern_x = -kern_pad; kern_x <= kern_pad; kern_x++) {       //run the kernel on one axis
				for (int kern_y = -kern_pad; kern_y <= kern_pad; kern_y++) {    //run the kernel on another axis
																				   //calculate the sum of points + kernel values surrounding the current point
					matrix_sum += src.at<uchar>(x + kern_x, y + kern_y)*sobel_horizontal.at<schar>(kern_x + kern_pad, kern_y + kern_pad);
				}
			}
			horizontal.at<int>(x, y) = (int)matrix_sum;
		}
	}

	Mat temp = Mat(src.rows, src.cols, CV_32S, Scalar(0));
	uint32_t max = 0;
	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			temp.at<int>(x, y) = (int32_t) sqrt(pow(vertical.at<int>(x, y), 2) + (pow(horizontal.at<int>(x, y), 2)));
			max = (max >= temp.at<int>(x, y) ? max : temp.at<int>(x, y));
			
		}
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));
	uint threshold = round(max*percent_threshold);
	
	for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
		for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

			if (temp.at<int>(x, y)>=threshold){
				out.at<uchar>(x, y) = 255;
			}
			else {
				out.at<uchar>(x, y) = 0;
			}
		}
	}

	//Mat angle = Mat(src.rows, src.cols, CV_16SC1, Scalar(0));
	//for (uint x = kern_pad; x < (src.rows - kern_pad); x++) {           //run the image on one axis
	//	for (uint y = kern_pad; y < (src.cols - kern_pad); y++) {      //run the image on another axis

	//		
	//			out.at<short>(x, y) =(int) round((atan(vertical.at<int>(x, y) / horizontal.at<int>(x, y))*180.0f)/ M_PI);

	//	}
	//}

	return out;
}

Mat vcpi_binary_blob_labelling(Mat src) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = src.clone();

	uint8_t label = 1;
	uint neigh_pixels[4];

	for (uint y = 1; y < src.rows - 1; y++) {
		for (uint x = 1; x < src.cols - 1; x++) {

			neigh_pixels[0] = out.at<uchar>(y - 1,x-1);
			neigh_pixels[1] = out.at<uchar>(y - 1, x);
			neigh_pixels[2] = out.at<uchar>(y - 1, x + 1);
			neigh_pixels[3] = out.at<uchar>(y, x - 1);

			if (src.at<uchar>(y, x) == 255)
			{
				if ((neigh_pixels[0] == 0) && (neigh_pixels[1] == 0 ) && (neigh_pixels[2]==0) && (neigh_pixels[3] == 0)) {
					out.at<uchar>(y, x) = label;
					label++;
				}
				else {
					uint8_t temp = 255;
					for (uint8_t i = 0; i < 4; i++)//finding the minimum element
					{
						if ((neigh_pixels[i] != 0) && (neigh_pixels[i] < temp)) {
							temp = neigh_pixels[i];
						}
					}
					out.at<uchar>(y, x) = temp;
					
					for (uint8_t i = 0; i < 4; i++)//replacing other elements
					{
						if ((neigh_pixels[i] < 255) && (neigh_pixels[i] > temp)) {
							out = find_replace_value(out, neigh_pixels[i],temp);
						}
					}
				}
			}
		}
	}
	return out;
}
Mat vcpi_binary_blob_improved_labelling(Mat src) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}

	Mat input = src.clone(); //input matrix to be analyzed
	Mat out = Mat(input.rows, input.cols, CV_8UC1, Scalar(0)); //output matrix
	Mat mask = Mat(input.rows, input.cols, CV_8UC1, Scalar(0)); //mask matrix
	uint8_t label = 1;

	for (uint y = 1; y < input.rows - 1; y++) {
		for (uint x = 1; x < input.cols - 1; x++) {

			if (input.at<uchar>(y, x) == 255) {
				mask = vcpi_gray_to_binary_Region_Growing(input,y,x, 10, 10, Oito);
				input ^= mask; //remove masked blob
				out |= find_replace_value(mask,255,label); //label the said blob
				label++; //increment labels
			}
		}
	}

	return out;
}

Mat vcpi_get_laser_line(Mat src) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}

	//remove the red channel to avoid artifacts
	//convert to grayscale
	Mat temp = vcpi_rgb_to_gray(vcpi_rgb_remove_red(src));

	temp = vcpi_gray_edge_sobel(temp, 0.7f);
	temp = vcpi_binary_close(temp, 20, Oito, 20);

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	uint up_pos = 0;
	uint low_pos = 0;

	for (uint x = 0; x < temp.cols; x++) { //find the edges of the line and get the middle position
		up_pos = 0;
		low_pos = 0;
		for (uint y = 0; y < temp.rows - 1; y++) {

			if ((temp.at<uchar>(y, x) == 0) && (temp.at<uchar>(y + 1, x) == 255)) {
				up_pos = y;
			}
			else if ((temp.at<uchar>(y, x) == 255) && (temp.at<uchar>(y + 1, x) == 0)) {
				low_pos = y;
			}
		}

		uint y_pos = (uint)ceil((up_pos + low_pos)*0.5f); //put the middle point of the line at 255
		out.at<uchar>(y_pos, x) = 255;
	}

	return out;
}

coordinates vcpi_blob_centroid(Mat src) {

	coordinates centroid = { centroid.x = 0, centroid.y = 0 };

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return centroid;
	}
	uint64_t area = 0;
	uint centerx = 0;
	uint centery = 0;

	for (uint y = 0; y < src.rows; y++) {
		for (uint x = 0; x < src.cols; x++) {

			if (src.at<uchar>(y, x) == 255) {

				centerx +=  x;
				centery += y;
				area++;
			}

		}
	}
		
	centroid.x = round(((double)centerx)/ ((double)area));
	centroid.y = round(((double)centery) / ((double)area));
	centroid.area = area;
	return centroid;
}

Mat vcpi_expanded_countour(Mat src, uint countour_thickness) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}

	Mat out = vcpi_binary_close(src, countour_thickness, Oito, 0, Quatro); //dilates the image by countour_thickness

	out ^= src;

	return out;
}

Mat vcpi_draw_circle_centroid(Mat src, uint circle_radius) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (circle_radius<=0)
	{
		cout << "Circle Radius must be positive!" << endl;
		return src;
	}
	
	coordinates centroid = vcpi_blob_centroid(src);

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	for (double y = 1; y < out.rows - 1; y++) {
		for (double x = 1; x < out.cols - 1; x++) {

			out.at<uchar>(y, x) = ((round(sqrt(pow((x - (double)centroid.x),2) + pow((y - (double)centroid.y),2))) <= circle_radius) ? 255 : 0);
		}
	}

	return out;
}


Mat vcpi_draw_line_between_centroids(Mat src, coordinates centroid_1, coordinates centroid_2, uint line_thickness) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (line_thickness <= 0)
	{
		cout << "Line Thickness must be positive!" << endl;
		return src;
	}

	uint half_thick_one = (uint)line_thickness / 2;
	uint half_thick_two = line_thickness % 2 == 0 ? half_thick_one : half_thick_one + 1;
	

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	double mx_normal = 0;
	double dividend_normal = (double)centroid_2.y - (double)centroid_1.y;
	double divider_normal  = round((double)centroid_2.x - (double)centroid_1.x);
	mx_normal = divider_normal == 0 ? 0 : (double)dividend_normal / (double)divider_normal;
	int b_normal = centroid_1.y - (centroid_1.x*mx_normal);

		double mx_inv = 0;
		double dividend_inv = (double)centroid_2.x - (double)centroid_1.x;
		double divider_inv = round((double)centroid_2.y - (double)centroid_1.y);
		mx_inv = divider_inv == 0 ? 0 : (double)dividend_inv / (double)divider_inv;
		int b_inv = centroid_1.x - (centroid_1.y*mx_inv);
	
	for (double y = 1; y < out.rows - 1; y++) {
		for (double x = 1; x < out.cols - 1; x++) {
			if (((x<=centroid_1.x && x >=centroid_2.x) || (x >= centroid_1.x && x <= centroid_2.x)) && ((y <= centroid_1.y && y >= centroid_2.y) || (y >= centroid_1.y && y <= centroid_2.y)))
			{
				if ((divider_normal != 0) && (abs(round(mx_normal)) <= 1))
				{
					double equation_normal = ((double)x * (double)mx_normal) + (double)b_normal;
					out.at<uchar>(y, x) = ((y >= (floor(equation_normal) - half_thick_one)) && (y <= (floor(equation_normal) + half_thick_two ))) ? 255 : 0;
				}
				else {
					double equation_inv = ((double)y * (double)mx_inv) + (double)b_inv;
					out.at<uchar>(y, x) = ((x >= (floor(equation_inv) - half_thick_two)) && (x <= (floor(equation_inv) + half_thick_two))) ? 255 : 0;
				}
				
			}
		}
	}

	return out;
}

coordinates home_point;

int point_compare(const void *vp1, const void *vp2)
{
	Point *p1 = (Point *)vp1;
	Point *p2 = (Point *)vp2;

	// Find orientation 
	int o;

	int val = ((*p1).y - home_point.y) * ((*p2).x - (*p1).x) - ((*p1).x - home_point.x) * ((*p2).y - (*p1).y);
	if (val == 0) {
		o = 0;
	}
	else {
		o = (val > 0) ? 1 : 2;
	}

	if (o == 0) {
		uint distance_one = (home_point.x - (*p1).x)*(home_point.x - (*p1).x) + (home_point.y - (*p1).y)*(home_point.y - (*p1).y);
		uint distance_two = (home_point.x - (*p2).x)*(home_point.x - (*p2).x) + (home_point.y - (*p2).y)*(home_point.y - (*p2).y);
		return (distance_two >= distance_one) ? -1 : 1;
	}

	return (o == 2) ? -1 : 1;
}

Mat vcpi_draw_line_labels_centroid(Mat src, uint line_thickness) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	else if (line_thickness <= 0)
	{
		cout << "Line Thickness must be positive!" << endl;
		return src;
	}

	uint max_labels = 0;
	for (uint i = 0; i < src.cols*src.rows; i++) {
		max_labels = src.data[i] > max_labels ? src.data[i] : max_labels;
	}

	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	if (max_labels>=2)
	{
	
	coordinates *points = (coordinates *)malloc((max_labels) * sizeof(coordinates));
	
	for (uint i = 0; i < max_labels; i++) {
		Mat blob = find_replace_value(src, i+1, 255);
		points[i] = vcpi_blob_centroid(blob);
	}

	int ymin = points[0].y, min = 0;
	for (uint i = 0; i < max_labels; i++)
	{
		// Pick the bottom-most. In case of tie, chose the 
		// left most point 
		if ((points[i].y < ymin) || (ymin == points[i].y && points[i].x < points[min].x)) {
			ymin = points[i].y;
			min = i;
		}
	}

	// Place the bottom-most point at first position 

	coordinates temp = points[0];
	points[0] = points[min];
	points[min] = temp;

	// Sort n-1 points with respect to the first point. 
	// A point p1 comes before p2 in sorted ouput if p2 
	// has larger polar angle (in counterclockwise 
	// direction) than p1 
	home_point = points[0];
	
	qsort(&points[1], max_labels - 1, sizeof(coordinates), point_compare);

	// Now stack has the output points, print contents 
	// of stack 

	
		for (uint i = 0; i < max_labels-1; i++)
		{
			out |= vcpi_draw_line_between_centroids(out, points[i], points[i + 1],line_thickness);
		}
		out |= vcpi_draw_line_between_centroids(out, points[0], points[max_labels-1], line_thickness);
	
	//draw line from the remaining point to the first point
	free(points);
	}
	return out;
}

Mat vcpi_filter_blob_area(Mat src,uint min_area) {

	if (src.empty()) {                	//check for input image
		cout << "There is no image!" << endl;
		return src;
	}
	Mat input = src.clone();
	Mat out = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));
	uint max_labels = 0;
	for (uint i = 0; i < src.cols * src.rows; i++) {
		max_labels = input.data[i] > max_labels ? input.data[i] : max_labels;
	}
	
	coordinates centroid;
	Mat temp = Mat(src.rows, src.cols, CV_8UC1, Scalar(0));

	for (uint i = 1; i <= max_labels; i++) {
		temp = find_replace_value(input, i, 255);
		centroid = vcpi_blob_centroid(temp);
		if (centroid.area >= min_area) {
			out |= temp;
		}
	}

	return out;
}