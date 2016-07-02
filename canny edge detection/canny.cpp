

#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>

#include<iostream>

///////////////////////////////////////////////////////////////////////////////////////////////////
int main() {

	cv::Mat imgOriginal;		// input image
	cv::Mat imgGrayscale;		// grayscale of input image in order to perform canny edge detection
	cv::Mat imgBlurred;			// intermediate blured image
	cv::Mat imgCanny;			// Canny edge image

	imgOriginal = cv::imread("image.jpg");			// open image

	if (imgOriginal.empty()) {									// if unable to open image
		std::cout << "error: image not read from file\n\n";		// show error message on command line
		return(0);												// and exit program
	}

	cv::cvtColor(imgOriginal, imgGrayscale, CV_BGR2GRAY);		// convert to grayscale

	cv::GaussianBlur(imgGrayscale,			// input image
		imgBlurred,							// output image
		cv::Size(7, 7),						// smoothing window width and height in pixels
		1.5);								// sigma value, determines how much the image will be blurred

	cv::Canny(imgBlurred,			// input image
		imgCanny,					// output image
		150,						// low threshold
		200);						// high threshold

									// declare windows
	cv::namedWindow("imgOriginal", CV_WINDOW_AUTOSIZE);	
	cv::namedWindow("imgCanny", CV_WINDOW_AUTOSIZE);		
															
	cv::imshow("imgOriginal", imgOriginal);		// show windows
	cv::imshow("imgCanny", imgCanny);

	cv::waitKey(0);					// hold windows open until user presses a key

	return(0);
}
