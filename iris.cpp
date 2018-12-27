#pragma once

#include "interface.h"
#include "iris.h"
#include "customize.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <cstdlib>

#define M_PI 3.1415926
static int _GetThresholdValudeFromOtsu(cv::Mat & matSrc);
using namespace cs;

int g_nStructElementsize = 3;

/**
 * 一种检测 pupilCircle 是否匹配瞳孔外圆的方法 
 */
static bool examine_pupilCircle(const Mat& img, const Circle& pupilCircle)
{
	if (
		pupilCircle.radius > img.cols / 5 ||
		pupilCircle.radius < img.cols / 15 ||
		get_distance(pupilCircle.center(), MyPoint(img.cols / 2, img.rows / 2)) > sqrt((pow(img.cols, 2) + pow(img.rows, 2))) / 5
		)
		return false;
	else
		return true;
}

static void _PrintPoint(cv::Mat &image, MyPoint &point)
{
	int xyadd[9][2] = {
		{ -1, -1 },
	{ -1, 0 },
	{ -1, 1 },
	{ 0, -1 },
	{ 0, 0 },
	{ 0, 1 },
	{ 1, -1 },
	{ 1, 0 },
	{ 1, 1 }
	};
	MyPoint _point = point;
	for (int i = 0; i < 9; i++)
	{
		_point = point;
		_point.x += xyadd[i][0];
		_point.y += xyadd[i][1];
		image.at<cv::Vec3b>(_point.y, _point.x)[0] = 0;
		image.at<cv::Vec3b>(_point.y, _point.x)[1] = 0;
		image.at<cv::Vec3b>(_point.y, _point.x)[2] = 255;
	}
}

//－－－－－－－－－－【用ＯＴＳＵ选择阀值】－－－－－－－－－－－－－－－－
static int _GetThresholdValudeFromOtsu(cv::Mat & matSrc)
{
	int nThreshold = 0;
	// 计算直方图 histogram
	float fHistogram[256] = { 0 };
	for (int y = 0; y < matSrc.rows; ++y)
	{
		uchar *pData = matSrc.ptr(y);
		for (int x = 0; x < matSrc.cols; ++x)
		{
			fHistogram[*pData++]++;
		}
	}
	int count = 0;
	for (int i = 200; i < 256; ++i)
	{
		count += fHistogram[i];
	}
	if (count < 200)
	{
		//flag = 1;
	}
	// 平均直方图 average histogram
	int nSize = matSrc.rows*matSrc.cols;
	for (int i = 0; i<256; i++)
	{
		fHistogram[i] = fHistogram[i] / nSize;
	}
	// 平均每个值 average pixel value
	float fAvaValue = 0;
	for (int i = 0; i<256; i++)
	{
		fAvaValue += i * fHistogram[i];
	}
	// 计算最大方差 max variance
	float fMaxVariance = 0;
	float fW = 0, fU = 0;
	for (int i = 0; i<256; i++)
	{
		fW += fHistogram[i];
		fU += i * fHistogram[i];
		float fTemp = fAvaValue * fW - fU;
		float fVariance = fTemp * fTemp / (fW*(1 - fW));
		if (fVariance > fMaxVariance)
		{
			fMaxVariance = fVariance;
			nThreshold = i;
		}
	}
	return nThreshold;
}


// return a valid pupilCircle(tough)
// 通过一个瞳孔内一点得到 粗原
static Circle get_tough_pupilCircle(const Mat& img, MyPoint point) {
	try {
		Mat dealingImg = img.clone();
		cv::GaussianBlur(dealingImg, dealingImg, Size(7, 7), 0);
		double minVal, maxVal;
		cv::minMaxLoc(dealingImg, &minVal, &maxVal);
		int minThreshold = minVal + 50;
		int maxThreshold = maxVal - 40;

		// liner finder
		auto _through_anchor_find_around =
			[&](const Mat& img, const MyPoint& point, int minThreshold, int maxThreshold) -> Around {
			Around around(point, point, point, point);

			int &min_x = around.left.x,
				&min_y = around.above.y,
				&max_x = around.right.x,
				&max_y = around.below.y;

			int cnt = 15;

			for (int val, count = 0, i = point.y; i > 0; i--)
				if (val = img.at<uchar>(i, point.x); val < minThreshold)
					min_y = i, count = 0;
				else if (!(count++ < cnt))
					break;

			for (int val, count = 0, i = point.y; i < img.rows; i++)
				if (val = img.at<uchar>(i, point.x); val < minThreshold)
					max_y = i, count = 0;
				else if (!(count++ < cnt))
					break;

			for (int val, count = 0, i = point.x; i > 0; i--)
				if (val = img.at<uchar>(point.y, i); val < minThreshold)
					min_x = i, count = 0;
				else if (!(count++ < cnt))
					break;

			for (int val, count = 0, i = point.x; i < img.cols; i++)
				if (val = img.at<uchar>(point.y, i); val < minThreshold)
					max_x = i, count = 0;
				else if (!(count++ < cnt))
					break;

			return around;
		};

		auto through_aroud_get_center = [](const Around& around) -> MyPoint {
			MyPoint center;
			center.x = (around.right.x + around.left.x) >> 1;
			center.y = (around.below.y + around.above.y) >> 1;
			return center;
		};


		Around around = _through_anchor_find_around(dealingImg, point, minThreshold, maxThreshold);

		MyPoint center = through_aroud_get_center(around);

		MyPoint reflectPoint(center.x + (center.x - point.x), (center.y + (center.y - point.y)));

		Around another_around = _through_anchor_find_around(dealingImg, reflectPoint, minThreshold, maxThreshold);

		vector<MyPoint> borders = through_around_get_points_vector({ around, another_around });

		//for (auto& point : borders) {
		//	draw_point(originalImage, point);
		//	draw_point(img, point);
		//}

		Circle _pupilCircle = circle_least_fit(borders);
		if (_pupilCircle.is_valid())
			return _pupilCircle;
		else
			return Circle(-1, -1, -1);
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		return Circle(-1, -1, -1);
	}
}


// 添加了另一个条件: val > 230
static Circle _get_tough_pupilCircle(const Mat& img, MyPoint point) {
	try {
		Mat dealingImg = img.clone();
		cv::GaussianBlur(dealingImg, dealingImg, Size(7, 7), 0);

		double minVal, maxVal;
		cv::minMaxLoc(dealingImg, &minVal, &maxVal);
		int minThreshold = minVal + 50;
		int maxThreshold = maxVal - 40;

		// liner finder
		auto _through_anchor_find_around =
			[&](const Mat& img, const MyPoint& point, int minThreshold, int maxThreshold) -> Around {
			Around around(point, point, point, point);

			int &min_x = around.left.x,
				&min_y = around.above.y,
				&max_x = around.right.x,
				&max_y = around.below.y;

			int cnt = 15;

			for (int val, count = 0, i = point.y; i > 0; i--)
				if (val = img.at<uchar>(i, point.x); val < minThreshold || val > 230)
					min_y = i, count = 0;
				else if (!(count++ < cnt))
					break;

			for (int val, count = 0, i = point.y; i < img.rows; i++)
				if (val = img.at<uchar>(i, point.x); val < minThreshold || val > 230)
					max_y = i, count = 0;
				else if (!(count++ < cnt))
					break;

			for (int val, count = 0, i = point.x; i > 0; i--)
				if (val = img.at<uchar>(point.y, i); val < minThreshold || val > 230)
					min_x = i, count = 0;
				else if (!(count++ < cnt))
					break;

			for (int val, count = 0, i = point.x; i < img.cols; i++)
				if (val = img.at<uchar>(point.y, i); val < minThreshold || val > 230)
					max_x = i, count = 0;
				else if (!(count++ < cnt))
					break;

			return around;
		};

		auto through_aroud_get_center = [](const Around& around) -> MyPoint {
			MyPoint center;
			center.x = (around.right.x + around.left.x) >> 1;
			center.y = (around.below.y + around.above.y) >> 1;
			return center;
		};


		Around around = _through_anchor_find_around(dealingImg, point, minThreshold, maxThreshold);

		MyPoint center = through_aroud_get_center(around);

		MyPoint reflectPoint(center.x + (center.x - point.x), (center.y + (center.y - point.y)));

		Around another_around = _through_anchor_find_around(dealingImg, reflectPoint, minThreshold, maxThreshold);

		vector<MyPoint> borders = through_around_get_points_vector({ around, another_around });

		//for (auto& point : borders) {
		//	draw_point(originalImage, point);
		//	draw_point(img, point);
		//}

		Circle _pupilCircle = circle_least_fit(borders);
		if (_pupilCircle.is_valid()) 
			return _pupilCircle;
		
		else
			return Circle(-1, -1, -1);
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		return Circle(-1, -1, -1);
	}
}


Iris::Iris(std::string imagePath) 
{
	this->imagePath = imagePath;
	add_image(imagePath);
}


void Iris::add_image(std::string imagePath) {
	this->imagePath = imagePath;
	try {

		grayImage = cv::imread(cv::String(imagePath), CV_LOAD_IMAGE_GRAYSCALE);
		originalImage = cv::imread(cv::String(imagePath));

		if (grayImage.empty() || originalImage.empty())
		{
			throw exception((string("can not open image at ") + imagePath).c_str());
		}
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
	}

}


MyPoint Iris::findPupilPoint()
{
	Mat img = grayImage.clone();
	MyPoint ret;

	//均值滤波
	blur(img, img, cv::Size(15, 15), cv::Point(-1, -1));

	//用Otsu对图形进行阀值分割
	unsigned char thre = _GetThresholdValudeFromOtsu(img);
	thre -= thre * 0.25;
	cv::threshold(img, img, thre, 255, cv::THRESH_BINARY);

	//对分割后的区域进行查找， 确定瞳孔内一点
	int k = 4;  //9*9区域进行查找
	int count_area_min = INT_MAX;
	for (int i = k; i < grayImage.rows - k; ++i)
	{
		for (int j = k; j < grayImage.cols - k; ++j)
		{
			if (img.at<uchar>(i, j) == 0)
			{
				//记录点附近灰度和
				int count_area = 0;
				for (int row = i - k; row <= i + k; ++row)
				{
					for (int col = j - k; col <= j + k; ++col)
					{
						count_area += grayImage.at<uchar>(row, col);
					}
				}
				//寻找灰度和最小的区域， 中心点作为要找的瞳孔区域内一点
				if (count_area < count_area_min)
				{
					count_area_min = count_area;
					ret = MyPoint(j, i);
				}
			}
		}
	}
	return ret;
}

// 原始的 test_way() 方法
bool Iris::original_find_pupilCircle()
{
	try
	{
		Mat img = grayImage.clone();
		//eliminate_border(img, img.cols / 6, img.rows / 6);
		MyPoint firstPoint = cs::find_pupil_inner_point(img);
		// 高斯滤波会干扰
		cv::GaussianBlur(img, img, Size(7, 7), 0);

		auto through_aroud_get_center = [](const Around& around) -> MyPoint {
			MyPoint center;
			center.x = (around.right.x + around.left.x) >> 1;
			center.y = (around.below.y + around.above.y) >> 1;
			return center;
		};

		double minVal, maxVal;
		cv::minMaxLoc(img, &minVal, &maxVal);

		int minThreshold = minVal + 50;
		int maxThreshold = maxVal - 40;

		// liner finder
		auto _through_anchor_find_around =
			[&](const Mat& img, const MyPoint& firstPoint, int minThreshold) -> Around {
			Around around(firstPoint, firstPoint, firstPoint, firstPoint);

			int &min_x = around.left.x,
				&min_y = around.above.y,
				&max_x = around.right.x,
				&max_y = around.below.y;


			for (int count = 0, i = firstPoint.y; i > 0; i--)
				if (img.at<uchar>(i, firstPoint.x) < minThreshold)
					min_y = i, count = 0;
				else if (!(count++ < 15))
					break;

			for (int count = 0, i = firstPoint.y; i < img.rows; i++)
				if (img.at<uchar>(i, firstPoint.x) < minThreshold)
					max_y = i, count = 0;
				else if (!(count++ < 15))
					break;

			for (int count = 0, i = firstPoint.x; i > 0; i--)
				if (img.at<uchar>(firstPoint.y, i) < minThreshold)
					min_x = i, count = 0;
				else if (!(count++ < 15))
					break;

			for (int count = 0, i = firstPoint.x; i < img.cols; i++)
				if (img.at<uchar>(firstPoint.y, i) < minThreshold)
					max_x = i, count = 0;
				else if (!(count++ < 15))
					break;

			return around;
		};

		// block finder: 待改进
		auto __through_anchor_find_around =
			[&](const Mat& img, const MyPoint& firstPoint, int minThreshold) -> Around {
			Around around(firstPoint, firstPoint, firstPoint, firstPoint);

			int &min_x = around.left.x,
				&min_y = around.above.y,
				&max_x = around.right.x,
				&max_y = around.below.y;

			auto pred = [&](uchar val) {return val < minThreshold; };
			int cnt = 100;
			int draft = 10;

			for (int count = 0, i = firstPoint.y; i > 0; i--) {

				MyPoint _leftUp(around.above.x - draft, around.above.y);
				MyPoint _rightDown(around.above.x + draft, around.above.y + draft);
				Mat_<uchar> roi = img(cv::Rect(_leftUp, _rightDown));

				if (std::count_if(roi.begin(), roi.end(), pred) > cnt
					) {
					min_y = i, count = 0;
					_leftUp.set(around.above.x - draft, around.above.y);
					_rightDown.set(around.above.x + draft, around.above.y + draft);
					roi = img(cv::Rect(_leftUp, _rightDown));

				}
				else if (!(count++ < 15))
					break;
			}

			for (int count = 0, i = firstPoint.y; i < img.rows; i++) {
				MyPoint _leftUp(around.below.x - draft, around.below.y - draft);
				MyPoint _rightDown(around.below.x + draft, around.below.y);
				Mat_<uchar> roi = img(cv::Rect(_leftUp, _rightDown));

				if (std::count_if(roi.begin(), roi.end(), pred) > cnt) {
					max_y = i, count = 0;
					_leftUp.set(around.below.x - draft, around.below.y - draft);
					_rightDown.set(around.below.x + draft, around.below.y);
					roi = img(cv::Rect(_leftUp, _rightDown));

				}
				else if (!(count++ < 15))
					break;
			}

			for (int count = 0, i = firstPoint.x; i > 0; i--) {
				MyPoint _leftUp(around.left.x, around.left.y - draft);
				MyPoint _rightDown(around.left.x + draft, around.left.y + draft);
				Mat_<uchar> roi = img(cv::Rect(_leftUp, _rightDown));

				if (std::count_if(roi.begin(), roi.end(), pred) > cnt) {
					min_x = i, count = 0;
					_leftUp.set(around.left.x, around.left.y - draft);
					_rightDown.set(around.left.x + draft, around.left.y + draft);
					roi = img(cv::Rect(_leftUp, _rightDown));

				}
				else if (!(count++ < 15))
					break;
			}

			for (int count = 0, i = firstPoint.x; i < img.cols; i++) {
				MyPoint _leftUp(around.left.x - draft, around.left.y - draft);
				MyPoint _rightDown(around.left.x, around.left.y + draft);
				Mat_<uchar> roi = img(cv::Rect(_leftUp, _rightDown));

				if (std::count_if(roi.begin(), roi.end(), pred) > cnt) {
					max_x = i, count = 0;
					_leftUp.set(around.left.x - draft, around.left.y - draft);
					_rightDown.set(around.left.x, around.left.y + draft);
					roi = img(cv::Rect(_leftUp, _rightDown));

				}
				else if (!(count++ < 15))
					break;
			}
			return around;
		};



		Around around = _through_anchor_find_around(img, firstPoint, minThreshold);

		MyPoint center = through_aroud_get_center(around);

		MyPoint reflectPoint(center.x + (center.x - firstPoint.x), (center.y + (center.y - firstPoint.y)));

		Around another_around = _through_anchor_find_around(img, reflectPoint, minThreshold);

		vector<MyPoint> borders = through_around_get_points_vector({ around, another_around });

		//for (auto& point : borders) {
		//	draw_point(originalImage, point);
		//	draw_point(img, point);
		//}

		pupilCircle = circle_least_fit(borders);

		if (pupilCircle.is_valid())
		{
			//if (examine_pupilCircle(grayImage, pupilCircle) == false) {
			//	cs::draw_circle(img, pupilCircle);
			//	cv::imshow(getImageName(), img);

			//	return false;
			//}
			//cv::circle(originalImage, pupilCircle.center(), pupilCircle.radius, Scalar(0, 0, 255), 2);
			//cv::imshow(getImageName(), originalImage);
			return true;
		}
		else
			return false;
	}
	catch (const exception& e)
	{
		cout <<"segment pupilCirlce :" << getImageName() << endl;
		cerr << e.what() << endl;
		return false;
	}
	return true;
}

// 你要的抠出来的东西
bool Iris::find_pupilCircle() {
	// common begin
	Mat dealingImg = grayImage.clone();
	eliminate_border(dealingImg, dealingImg.cols / 6, dealingImg.rows / 6);
	// common end

	MyPoint firstPoint = cs::find_pupil_inner_point(dealingImg);
	
	// 判断得到的粗圆 是否正确: 通过 .is_valid 方法判断
	pupilCircle = get_tough_pupilCircle(dealingImg, firstPoint);

	//MyPoint secondPoint = findPupilPoint();
	//pupilCircle = get_tough_pupilCircle(dealingImg, secondPoint);

	if (pupilCircle.is_valid())
		return true;
	else
		return false;
}

bool Iris::getPupilCircle() noexcept 
{
	return findPupilEdge();
}

pair<MyPoint, MyPoint> Iris::getPupilCenterArea(int xRadius, int yRadius) {
	if (find_pupilCircle()) {
		return make_pair(
			MyPoint(pupilCircle.centerX - xRadius, pupilCircle.centerY - yRadius),
			MyPoint(pupilCircle.centerX + xRadius, pupilCircle.centerY + yRadius)
		);
	}
	else {
		throw exception("error occur since pupilCircle is invalid");
	}
}

//计算出瞳孔中心点及半径
bool Iris::findPupilEdge()
{
	int px = 0, py = 0, pr = 0;

	int range_x = 8, range_y = 8;
	pair<MyPoint, MyPoint> 	 twoPoint;
	try {
		twoPoint = getPupilCenterArea(range_x, range_y);
	}
	catch (exception& e) {
		cerr << e.what() << endl;
		return false;
	}

	MyPoint leftAbove = twoPoint.first;
	MyPoint rightBelow = twoPoint.second;
	MyPoint center = MyPoint(
		leftAbove.x + 8,
		leftAbove.y + 8
	);
	int width = 16, height = 16;

	//int height = range_x ;
	// int width = 32;
	//int ratio = 1;
	int max = 0;
	//range_x *= ratio;
	//range_y *= ratio;
	for (int i = leftAbove.y; i < leftAbove.y + width; ++i)
		// for (int i = center.y - range_y; i < center.y + range_y; ++i)
	{
		//for(int j = leftAbove.x;j < leftAbove.x+height; ++j)
		for (int j = center.x - range_x; j < center.x + range_x; ++j)
		{
			for (int r = 20; r < 60; r++)
			{
				int die = calCircleSum_pupil(grayImage, j, i, r + 1) - calCircleSum_pupil(grayImage, j, i, r - 1);
				if (die > max)
				{
					px = j; py = i; pr = r;
					max = die;
				}
			}
		}
	}

	//for (int i = leftAbove.y; i < leftAbove.y + width; ++i)
	//{
	//	for (int j = leftAbove.x; j < leftAbove.x + height; ++j)
	//	{
	//		grayImage.at<uchar>(i, j) = 0;
	//	}
	//}
	//cv::imshow(getImageName() + "0", grayImage);
	pupilCircle.centerX = px;
	pupilCircle.centerY = py;
	pupilCircle.radius = pr;

	cv::circle(originalImage, pupilCircle.center(), pupilCircle.radius, Scalar(0, 0, 255), 2);

	return true;
}


int Iris::calCircleSum_pupil(cv::Mat& img, int x, int y, int r) {
	auto radian = [=](double degree) {
		static const double PI = 3.1415926;
		return degree / (180 / PI);
	};
	unsigned sum = 0;
	double alpha;
	int cx, cy;
	pair<int, int> range{ -180, 180 };
	for (int degree = range.first; degree <= range.second; ++degree)
	{

		alpha = radian(degree);
		cx = x + r * cos(alpha);
		cy = y + r * sin(alpha);
		// 使得 cy, cx 作为索引不超出图像所能索引的最大范围
		// cy 
		if (cy > img.rows - 1)
			cy = img.rows - 1;
		else if (cy < 0)
			cy = 0;
		uchar *data = img.ptr<uchar>(cy);
#if defined(VERBOSE)
		//  cx
		if (cx + 1 > img.cols - 1)
			cx = img.cols - 2;
		else if (cx - 1 < 0)
			cx = 1;
		sum += (unsigned)data[cx - 1] + (unsigned)data[cx] + (unsigned)data[cx + 1];
#else	
		if (cx > img.cols - 1)
			cx = img.cols - 1;
		if (cx < 0)
			cx = 0;
		sum += data[cx];
#endif // VERBOSE
	}
	return (int)sum;

}

int Iris::calCircleSum(cv::Mat &img, int x, int y, int r)
{
#define VERBOSE // 用三圈
#define FORMUAL_EXPR
#if defined(FORMUAL_EXPR) // 公式中的方法
	auto radian = [=](double degree) {
		static const double PI = 3.1415926;
		return degree / (180 / PI);
	};
	unsigned sum = 0, sums[2] = { 0, 0 };
	double alpha;
	int cx, cy;
	vector<pair<int, int>> degree_ranges{ {-45, 45}, {135, 225}, {-30, 30}, {160, 210} };

	for (int s = 0; s < 2; ++s)
	for (int i = s * 2; i < 2 + s * 2; ++i)
		for (const auto & range : degree_ranges)
			for (int degree = range.first; degree <= range.second; ++degree)
			{

				alpha = radian(degree);
				cx = x + r * cos(alpha);
				cy = y + r * sin(alpha);
				// 使得 cy, cx 作为索引不超出图像所能索引的最大范围
				// cy 
				if (cy > img.rows - 1)
					cy = img.rows - 1;
				else if (cy < 0)
					cy = 0;
				uchar *data = img.ptr<uchar>(cy);
#if defined(VERBOSE)
				//  cx
				if (cx + 1 > img.cols - 1)
					cx = img.cols - 2;
				else if (cx - 1 < 0)
					cx = 1;
				sums[s%2] += (unsigned)data[cx - 1] + (unsigned)data[cx] + (unsigned)data[cx + 1];
#else	
				if (cx > img.cols - 1)
					cx = img.cols - 1;
				if (cx < 0)
					cx = 0;
				sums[s % 2] += (unsigned)data[cx];
#endif
			}
	sum = (sums[0] > sums[1]) ? sums[0]: sums[1];
#else // 刚才用的
	auto radian = [=](double degree) {
		static const double PI = 3.1415926;
		return degree / (180 / PI);
	};
	unsigned sum = 0;
	double alpha;
	int cx, cy;
	vector<pair<int, int>> degree_ranges{ {-45, 45}, {135, 225}, {-30, 30}, {160, 210} };

	for (int i = 0; i < 2 + degree_ranges.size(); ++i)
		for (const auto & range : degree_ranges)
			for (int degree = range.first; degree <= range.second; ++degree)
			{

				alpha = radian(degree);
				cx = x + r * cos(alpha);
				cy = y + r * sin(alpha);
				// 使得 cy, cx 作为索引不超出图像所能索引的最大范围
				// cy 
				if (cy > img.rows - 1)
					cy = img.rows - 1;
				else if (cy < 0)
					cy = 0;
				uchar *data = img.ptr<uchar>(cy);
#if defined(VERBOSE)
				//  cx
				if (cx + 1 > img.cols - 1)
					cx = img.cols - 2;
				else if (cx - 1 < 0)
					cx = 1;
				sum += (unsigned)data[cx - 1] + (unsigned)data[cx] + (unsigned)data[cx + 1];
#else	
				if (cx > img.cols - 1)
					cx = img.cols - 1;
				if (cx < 0)
					cx = 0;
				sums += data[cx];
#endif // VERBOSE
			}
#endif // FORMUAL_EXPR
	return (int)sum;
}

/*
bool Iris::findIrisEdge(const cv::Mat &src)
{
	try {
		int imax = 0;
		int iris_r;

		cv::GaussianBlur(src, src, Size(5, 5), 0);
		for (int r1 = pupilCircle.radius * 1.5; r1 < pupilCircle.radius * 3.5; r1++)
		{
			int die1 = calCircleSum(src, pupilCircle.centerX, pupilCircle.centerY, r1 + 1) - calCircleSum(src, pupilCircle.centerX, pupilCircle.centerY, r1 - 1);
			if (die1 > imax) {
				iris_r = r1;
				imax = die1;
			}
		}
		irisCircle = Circle(pupilCircle.centerX, pupilCircle.centerY, iris_r - 1);
	}
	catch (const exception& e) {
		return false;
	}
	if (irisCircle.is_valid())
		return true;
	else
		return false;
}
*/
/*
// @parm src: 原图的灰度图像
bool Iris::findIrisEdge(const cv::Mat &_src)
{
	try {
		Mat src = _src.clone();
		cv::GaussianBlur(src, src, Size(5, 5), 0);
		int imax = 0;
		int iris_r;
		// semi_slide : pupilCircle 内接矩形的边长的一半
		//int semi_slide = double(sqrt(2) * .25 * double(pupilCircle.radius)); // pupulCircle 已经找到的瞳孔外圆(虹膜内圆)
		int semi_slide = 8;

		// 对内接正方形的所有点进行遍历
		for (int i = pupilCircle.center().x - semi_slide; i < pupilCircle.center().x + semi_slide; ++i)
		{
			for (int j = pupilCircle.center().y - semi_slide; j < pupilCircle.center().y + semi_slide; ++j)
			{
				for (int r1 = pupilCircle.radius * 1.5; r1 < pupilCircle.radius * 3; r1++)
				{
					//int die1 = calCircleSum(src, pupilCircle.centerX, pupilCircle.centerY, r1 + 1) - calCircleSum(src, pupilCircle.centerX, pupilCircle.centerY, r1 - 1);
					int die1 = calCircleSum(src, i, j, r1 + 1) - calCircleSum(src, i, j, r1 - 1);
					if (die1 > imax) {
						iris_r = r1;
						imax = die1;
						irisCircle.centerX = i;
						irisCircle.centerY = j;
					}
				}
			}
		}
		int _;
		irisCircle = Circle(irisCircle.centerX, irisCircle.centerY, iris_r - 1); // irisCircle 是虹膜外圆
	}
	catch (const exception& e) {
		return false;
	}
	// if(is_valid == true) --> isrCircle.centerX > 0 && irisCircle.centerY > 0 && irisCircle.radius
	// else --> false
	if (irisCircle.is_valid())
		return true;
	else
		return false;
} // return false 找圆失败
*/

#ifdef SAME_CENTER
// 9-29 17:44
bool Iris::findIrisEdge(const cv::Mat &_src)
{
	try {
		Mat src = _src.clone();
		cv::GaussianBlur(src, src, Size(7, 7), 1.0);

		int g_nStructElementsize = 3;
		cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE,
			cv::Size(2 * g_nStructElementsize + 1, 2 * g_nStructElementsize + 1),
			cv::Point(g_nStructElementsize, g_nStructElementsize));
		morphologyEx(src, src, cv::MORPH_CLOSE, element);
		morphologyEx(src, src, cv::MORPH_OPEN, element);



		int imax = 0;
		int iris_r;
		// semi_slide : pupilCircle 内接矩形的边长的一半
		int semi_slide = double(sqrt(2) * .25 * double(pupilCircle.radius)); // pupulCircle 已经找到的瞳孔外圆(虹膜内圆)

		// 对内接正方形的所有点进行遍历

		for (int i = pupilCircle.radius - semi_slide; i < pupilCircle.radius + semi_slide; ++i)
		for (int r1 = pupilCircle.radius * 1.5; r1 < pupilCircle.radius * 3.5; r1++)
		{
			long long die1 = (
				calCircleSum(src, pupilCircle.centerX, pupilCircle.centerY, r1 + 1)
				- calCircleSum(src, pupilCircle.centerX, pupilCircle.centerY, r1 - 1)
				);
			//int die1 = calCircleSum(src, i, j, r1 + 1) - calCircleSum(src, i, j, r1 - 1);
			if (die1 > imax) {
				iris_r = r1;
				imax = die1;
			}
		}
		irisCircle = Circle(pupilCircle.centerX, pupilCircle.centerY, iris_r - 1); // pupilCircle 是虹膜外圆
	}
	catch (const exception& e) {
		return false;
	}
	// if(is_valid == true) --> isrCircle.centerX > 0 && irisCircle.centerY > 0 && irisCircle.radius
	// else --> false
	if (irisCircle.is_valid())
		return true;
	else
		return false;
} // return false 找圆失败
#elif defined(DIFF_CENTER)

bool Iris::findIrisEdge(const cv::Mat &_src)
{
	try {
		Mat src = _src.clone();
		cv::GaussianBlur(src, src, Size(7, 7), 1.0);

		int g_nStructElementsize = 3;
		cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE,
			cv::Size(2 * g_nStructElementsize + 1, 2 * g_nStructElementsize + 1),
			cv::Point(g_nStructElementsize, g_nStructElementsize));

		morphologyEx(src, src, cv::MORPH_CLOSE, element);
		morphologyEx(src, src, cv::MORPH_OPEN, element);

		int imax = 0;
		int iris_r;
		int irisX, irisY;
		// semi_slide : pupilCircle 内接矩形的边长的一半
		int semi_slide = 5; // double(sqrt(2) * .25 * double(pupilCircle.radius)); // pupulCircle 已经找到的瞳孔外圆(虹膜内圆)
			 
		// 对内接正方形的所有点进行遍历

		for (int i = pupilCircle.centerX - semi_slide; i < pupilCircle.centerX + semi_slide; ++i)
			for (int j = pupilCircle.centerY - semi_slide; j < pupilCircle.centerY + semi_slide; ++j)
				for (int r1 = pupilCircle.radius * 1.5; r1 < pupilCircle.radius * 4; r1++)
				{
					long long die1 = (
						calCircleSum(src, i, j, r1 + 1)
						- calCircleSum(src, i, j, r1 - 1)
						);
					if (die1 > imax) {
						iris_r = r1;
						imax = die1;
						irisX = i;
						irisY = j;
					}
				}
		irisCircle = Circle(irisX, irisY, iris_r - 1); // pupilCircle 是虹膜外圆
	}
	catch (const exception& e) {
		return false;
	}
	// if(is_valid == true) --> isrCircle.centerX > 0 && irisCircle.centerY > 0 && irisCircle.radius
	// else --> false
	if (irisCircle.is_valid())
		return true;
	else
		return false;
} // return false 找圆失败
#endif

bool Iris::getIrisCircle() noexcept
{
	if (findIrisEdge(grayImage))
		return true;
	else
		return false;
}

void Iris::drawCircle() noexcept {
	if (pupilCircle.is_valid())
		cv::circle(originalImage, pupilCircle.center(), pupilCircle.radius, Scalar(0, 0, 255), 2);
	if (irisCircle.is_valid())
		cv::circle(originalImage, irisCircle.center(), irisCircle.radius, Scalar(0, 255, 255), 2);
	cv::imshow(getImageName(), originalImage);
}
