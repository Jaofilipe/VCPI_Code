#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;

typedef enum TipoQuadrado {  // typedef to set some names to help with shading
	left = 1,
	top = 2,
	corner = 3
} TipoQuadrado;

typedef enum TipoVizinhanca {  // typedef to set some names to help with shading
	Quatro = 4,
	Oito = 8
} TipoVizinhanca;

extern Mat drawShadedSquare(TipoQuadrado tipo);
extern Mat find_replace_value(Mat src, uint value, uint replace_with);
extern Mat vcpi_gray_negative(Mat src);
extern Mat vcpi_rgb_negative(Mat src);
extern Mat vcpi_rgb_remove_red(Mat src);
extern Mat vcpi_rgb_remove_blue(Mat src);
extern Mat vcpi_rgb_remove_green(Mat src);
extern Mat vcpi_rgb_to_gray(Mat src);
extern Mat vcpi_rgb_get_R(Mat src);
extern Mat vcpi_rgb_get_G(Mat src);
extern Mat vcpi_rgb_get_B(Mat src);
extern Mat vcpi_gray_to_binary(Mat src, int Lower_threshold = 0, int Upper_threshold = 255);
extern Mat VCPI_Segmenta_Cor(Mat src, int Lower_h = 0, int Upper_h = 255, int Lower_s = 0,
	int Upper_s = 255, int Lower_v = 0, int Upper_v = 255);
extern Mat vcpi_scale_gray_to_rgb(Mat src);
extern Mat vcpi_convolucao(Mat src, Mat kernel);
extern Mat vcpi_median_filter(Mat src, uint kernel_size = 3);
extern Mat vcpi_gray_to_binary_global_mean(Mat src);
extern Mat vcpi_gray_to_binary_midpoint(Mat src, uint kernel_size = 3);
extern Mat vcpi_gray_to_binary_bernsen(Mat src, uint kernel_size = 3, uint Cmin = 3);
extern Mat vcpi_gray_to_binary_niblack(Mat src, uint kernel_size = 3, float k = 1.0f);
extern Mat vcpi_gray_to_binary_Region_Growing(Mat src, uint x, uint y, uint Lower_Limit = 10,
	uint Upper_Limit = 10, TipoVizinhanca neighborhood = Quatro);
extern Mat vcpi_gray_to_binary_otsu(Mat src);