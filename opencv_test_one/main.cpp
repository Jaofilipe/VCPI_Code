
#include <opencv2/highgui.hpp>
#include <iostream>
using namespace cv;
using namespace std;
int main()
{
	Mat im = imread("C:/Users/João Filipe/Pictures/jets.png");
	if (im.empty()) {
		cout << "Cannot load file" << endl;
		return -1;
	}

	cout << "Channels:" << im.channels() << endl;
	cout << "Rows:" << im.rows << endl;
	cout << "Columns:" << im.cols << endl;

	for (int i = 0; i < (im.rows*im.cols*im.channels());i+=im.channels()) {
		im.data[i] = 255 - im.data[i];
		im.data[i + 1] = im.data[i + 1];
		im.data[i + 2] = 255 - im.data[i + 2];
	}

	imshow("Image1", im);

	waitKey(0);

	return 0;
}