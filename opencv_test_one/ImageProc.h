#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"

#define M_PI 3.14159265358979323846

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

typedef struct coordenadas {
	uint64_t x;
	uint64_t y;
	uint64_t area;
}coordinates;

extern Mat drawShadedSquare(TipoQuadrado tipo);
extern Mat find_replace_value(Mat src, uint value, uint replace_with);
extern Mat find_create_with(Mat src, uint value, uint red, uint green, uint blue);
extern Mat find_apply_created(Mat src, Mat src_created);
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
extern Mat vcpi_convolucao_freqs(Mat src, Mat kernel, double divide_by = 1);
extern Mat vcpi_median_filter(Mat src, uint kernel_size = 3);
extern Mat vcpi_gray_to_binary_global_mean(Mat src);
extern Mat vcpi_gray_to_binary_midpoint(Mat src, uint kernel_size = 3);
extern Mat vcpi_gray_to_binary_bernsen(Mat src, uint kernel_size = 3, uint Cmin = 3);
extern Mat vcpi_gray_to_binary_niblack(Mat src, uint kernel_size = 3, float k = 1.0f);
extern Mat vcpi_gray_to_binary_Region_Growing(Mat src, uint x, uint y, uint Lower_Limit = 10,
	uint Upper_Limit = 10, TipoVizinhanca neighborhood = Quatro);
extern Mat vcpi_gray_to_binary_otsu(Mat src);
extern Mat vcpi_dilation(Mat src, TipoVizinhanca neighborhood = Quatro);
extern Mat vcpi_erosion(Mat src, TipoVizinhanca neighborhood = Quatro);
extern Mat vcpi_binary_close(Mat src, uint dilation_times = 2, TipoVizinhanca dilation_neighborhood = Quatro, uint erosion_times = 2, TipoVizinhanca erosion_neighborhood = Quatro);
extern Mat vcpi_binary_open(Mat src, uint erosion_times = 2, TipoVizinhanca erosion_neighborhood = Quatro, uint dilation_times = 2, TipoVizinhanca dilation_neighborhood = Quatro);
extern Mat vcpi_gray_edge_prewitt(Mat src, float th = 0.5f);
extern Mat vcpi_gray_edge_sobel(Mat src, float th=0.5f);
extern Mat vcpi_binary_blob_labelling(Mat src);
extern Mat vcpi_binary_blob_improved_labelling(Mat src);
extern Mat vcpi_get_laser_line(Mat src);
extern coordinates vcpi_blob_centroid(Mat src);
extern Mat vcpi_expanded_countour(Mat src, uint countour_thickness = 1);
extern Mat vcpi_draw_circle_centroid(Mat src,uint circle_radius=4);
extern Mat vcpi_draw_line_between_centroids(Mat src, coordinates centroid_1, coordinates centroid_2, uint line_thickness = 2);
extern Mat vcpi_draw_line_labels_centroid(Mat src,uint line_thickness=2);
extern Mat vcpi_filter_blob_area(Mat src, uint min_area = 1);