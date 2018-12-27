#pragma once

#include "assistant.h"
#include "interface.h"
#include <initializer_list>


// 自定义的函数库
namespace customize {
	// 返回 Point(x1, y1) 和 Point(x2, y2) 之间的距离
	inline int get_distance(double x1, double y1, double x2, double y2)
	{
		return cvRound(sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2)));
	}

	// 返回 a 和 b 之间的距离
	inline int get_distance(const Point& a, const Point& b)
	{
		return get_distance(a.x, a.y, b.x, b.y);
	}

	/*
	store four points
	存储四个MyPoint对象, 这四个MyPoint对象应为
	瞳孔内 某一点 对应的 上下左右 四个四个瞳孔的边界点
	*/
	struct Around {
		MyPoint above, below, left, right;
		Around() = default;
		Around(const MyPoint& above, const MyPoint& below, const MyPoint& left, const MyPoint& right)
			:above(above), below(below), left(left), right(right)
		{}
		Around(const Around& around)
		{
			this->operator=(around);
		}

		void set_around(const MyPoint& above, const MyPoint& below, const MyPoint& left, const MyPoint& right)
		{
			this->above = above, this->below = below, this->left = left, this->right = right;
		}

		const Around& operator=(const Around& around)
		{
			above = around.above,
				below = around.below,
				left = around.left,
				right = around.right;
			return *this;
		}

		bool operator<(Around& around)
		{
			return above.y + left.x < around.above.y + around.left.x;
		}
	};

#ifdef BASE_DIR

	/*
	返回图片的绝对路径, 要求BASE_DIR已经定义过
	imgName 格式: "^S1\d{3}[LR]\d{2}.jpg$" -示例-> "S1023R02.jpg"
	get_abspath function have one parameter of imgName with type of std::string
	which will return the absolute path of the image in the computer
	notice that this fuction will use the macro of BASE_DIR which represent the root path of the images, thus BASE_DIR must defined beofore.
	*/
	std::string get_abspath(std::string imgName) {
		return std::string(BASE_DIR) + "//" + imgName.substr(2, 3) + "//" + imgName.substr(5, 1) + "//" + imgName;
	}
#endif

	// 返回图片的绝对路径
	// imgName 格式: "^S1\d{3}[LR]\d{2}.jpg$" -示例-> "S1023R02.jpg"
	std::string get_abspath(std::string imgName, std::string base_dir) {
		return base_dir + "//" + imgName.substr(2, 3) + "//" + imgName.substr(5, 1) + "//" + imgName;
	}

	/*
	计算区域内像素点的加和， 区域 = { Point(x, y) | startX <= x < col && startY <= y < row }
	accumlate fuction will return the accumlation of the img from the area of { (x, y) | startX <= x < col && startY <= y < row }
	*/
	int accumulate(const Mat& img, int startX, int startY, int col, int row)
	{
		int sum = 0;
		int endX = startX + col;
		int endY = startY + row;

		for (int i = startX; i < endX; ++i)
			for (int j = startY; j < endY; ++j)
				sum += img.at<uchar>(j, i);

		return sum;
	}

	/*
	计算区域内像素点的加和, 区域 = { Point(x, y) | startX <= x < col && startY <= y < row }
	计算 以startPoint为起点col列，row行区域内像素点的累加和
	*/
	int accumulate(const Mat& img, const Point& startPoint, int col, int row)
	{
		return accumulate(img, startPoint.x, startPoint.y, col, row);
	}

	/*
	用 value 填充 { Point(x, y) | startX <= x < col && startY <= y < row }
	fill_value fuction will assign value in every point at the area of { Point(x, y) | startX <= x < col && startY <= y < row }
	*/
	void fill_value(Mat& img, int value, int startX, int startY, int col, int row)
	{
		int endX = startX + col;
		int endY = startY + row;
		for (int i = startX; i < endX; ++i)
		{
			for (int j = startY; j < endY; ++j)
				img.at<uchar>(j, i) = cv::saturate_cast<uchar>(value);
		}
	}

	/*
	返回一个 Mat_<char> 类型的对象，可以被用来作为cv::filter2D函数的kernel参数
	可以有效的去除睫毛
	return a kernel which can be used in cv::filter2D as a kernel argument
	which aim for eliminate eyelash
	*/
	Mat_<char> centralizeKernel(int size = 5) {
		Mat_<char> kernel(Size(size, size));
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j) {
				kernel.at<char>(i, j) = (i == size / 2 || j == size / 2) ? 0 : -1;
			}
		kernel.at<char>(size / 2, size / 2) = size * size - 4;
		return move(kernel);
	}

	/*
	将图像img的左右xlen列, 上下ylen行 的值填充为val,
	val 默认为255 即白色，
	如果 ylen = 0,则ylen被赋值为 xlen
	*/
	void eliminate_border(Mat& img, int xlen, int ylen = 0, int val = 255)
	{
		if (ylen == 0) ylen = xlen;

		assert(xlen < img.cols / 2 && ylen < img.rows / 2);

		for (int i = 0; i < xlen; ++i)
			for (int j = 0; j < img.rows; ++j)
				img.at<uchar>(j, i) = val;
		for (int i = img.cols - 1; i >= img.cols - xlen; --i)
			for (int j = 0; j < img.rows; ++j)
				img.at<uchar>(j, i) = val;
		for (int i = 0; i < img.cols; ++i)
			for (int j = img.rows - 1; j >= img.rows - ylen; --j)
				img.at<uchar>(j, i) = val;
		for (int i = 0; i < img.cols; ++i)
			for (int j = 0; j < ylen; ++j)
				img.at<uchar>(j, i) = val;
	}

	/*
		将参数链表中的所有around中的所有四个点 放到一个vector中
	*/
	vector<MyPoint> through_around_get_points_vector(initializer_list<Around> arounds)
	{
		vector<MyPoint> pointsVec;
		for (Around around : arounds)
		{
			pointsVec.push_back(around.above);
			pointsVec.push_back(around.below);
			pointsVec.push_back(around.left);
			pointsVec.push_back(around.right);
		}
		return move(pointsVec);
	}


	/*
	先对grayImg进行拷贝为pyrImg
	先对pyrImg进行两次pyrDown操作，使得图像的的尺寸减小为原来的1/4
	再进行两次pyrUp操作 使得图像的尺寸恢复回去
	对变换后的图像(pyrImg)进行霍夫圆变换，返回变换后的 圆向量
	*/
	vector<cv::Vec3f> pixel_blur2(const Mat& grayImg)
	{
		Mat pyrImg;
		pyrDown(grayImg, pyrImg, Size(grayImg.cols / 2, grayImg.rows / 2));
		pyrDown(pyrImg, pyrImg, Size(pyrImg.cols / 2, pyrImg.rows / 2));

		pyrUp(pyrImg, pyrImg, Size(pyrImg.cols * 2, pyrImg.rows * 2));
		pyrUp(pyrImg, pyrImg, Size(pyrImg.cols * 2, pyrImg.rows * 2));

		vector<cv::Vec3f> circles;
		HoughCircles(pyrImg, circles, cv::HOUGH_GRADIENT, 2.4, 10, 100, 100);
		return circles;

	}

	/*
	先对grayImg进行拷贝为pyrImg
	先对pyrImg进行一次pyrDown操作，使得图像的的尺寸减小为原来的1/2
	再进行一次pyrUp操作 使得图像的尺寸恢复回去
	对变换后的图像(pyrImg)进行霍夫圆变换，返回变换后的 圆向量
	*/
	vector<cv::Vec3f> pixelBlur(const Mat& grayImg)
	{
		Mat pyrImg;
		pyrDown(grayImg, pyrImg, Size(grayImg.cols / 2, grayImg.rows / 2));
		pyrUp(pyrImg, pyrImg, Size(pyrImg.cols * 2, pyrImg.rows * 2));

		vector<cv::Vec3f> circles;
		HoughCircles(pyrImg, circles, cv::HOUGH_GRADIENT, 2.4, 10, 100, 100);
		return circles;
	}

	/*
	在img上用+号标记出point的位置，
	加号的半径可选，默认为3
	加号的颜色可选，默认为紫红色
	*/
	void draw_point(Mat& img, const Point& point, int radius = 3, Scalar color = Scalar(255, 0, 255))
	{
		cv::line(img, Point(point.x - radius, point.y), Point(point.x + radius, point.y), color);
		cv::line(img, Point(point.x, point.y - radius), Point(point.x, point.y + radius), color);
	}


	// 返回在srcImg 中 circle 圆内，像素值低于threshold的百分比
	double get_below_threshold_rate(const Mat& srcImg, const Circle& circle, int threshold = 100) {
		// double rate;
		int count = 0;
		int belowThresholdCount = 0, aboveThresholdCount = 0;
		for (int i = circle.centerY - circle.radius; i < circle.centerY + circle.radius; ++i)
		{
			for (int j = circle.centerX - circle.radius; j < circle.centerX + circle.radius; ++j)
			{
				if (circle.in_circle(i, j))
				{
					++count;
					if (srcImg.at<uchar>(i, j) < threshold)
						++belowThresholdCount;
				}

			}
		}
		return ((double)belowThresholdCount) / ((double)count);
	};

	/*
	对primitiveCircles中的圆进行检测，将符合特征的圆放入validCircles 容器内
	for each cirlce in primitiveCircles find circles which suit for the characteristics of pupil cirlce
	and put these circle to dealing vector
	*/
	void circle_screen(const Mat& srcImg, vector<cv::Vec3f>& primitiveCircles, vector<cv::Vec3f>& validCircles)
	{
		// 检测圆内 像素值小于 threshold 的点在 圆内占的比例
		auto innerCircleExamine = [&](Circle&& circle, int threshold = 100) {
			// double rate;
			int count = 0;
			int belowThresholdCount = 0, aboveThresholdCount = 0;
			for (int i = circle.centerY - circle.radius; i < circle.centerY + circle.radius; ++i)
			{
				for (int j = circle.centerX - circle.radius; j < circle.centerX + circle.radius; ++j)
				{
					if (circle.in_circle(i, j))
					{
						++count;
						if (srcImg.at<uchar>(i, j) < threshold)
							++belowThresholdCount;
					}

				}
			}
			return ((double)belowThresholdCount) / ((double)count);
		};

		MyPoint imgCenter(srcImg.cols / 2, srcImg.rows / 2);
		if (validCircles.size() != 0)
			validCircles.clear();

		for (int i = 0; i < primitiveCircles.size(); i++) {

			int centerX = cvRound(primitiveCircles.at(i)[0]);
			int centerY = cvRound(primitiveCircles.at(i)[1]);
			int radius = cvRound(primitiveCircles.at(i)[2]);
			MyPoint circleCenter = MyPoint(centerX, centerY);
			if (
				abs(circleCenter.x - imgCenter.x) < srcImg.cols / 4
				&& abs(circleCenter.y - imgCenter.y) < srcImg.rows / 4
				&& radius < srcImg.cols / 5
				&& radius > sqrt(pow(srcImg.cols, 2) + pow(srcImg.rows, 2)) / 18
				&& innerCircleExamine(Circle(primitiveCircles.at(i))) > .6
				)
			{
				validCircles.push_back(primitiveCircles.at(i));
			}
		}

	}

	/*
	和Mat::coverTo(Mat&) 的最用一样，调节图像的对比度
	same effect as the fucntion as Mat::convertTo
	*/
	void convert(Mat& src, Mat& dst, double alpha, double beta)
	{
		Mat_<uchar> table(Size(256, 1), CV_8UC1);
		for (int i = 0; i < 256; ++i)
			*(table.data + i) = ((i*alpha + beta) > 255 ? 255 : (i*alpha + beta) < 0 ? 0 : (i*alpha + beta));
		cv::LUT(src, table, dst);
	}

	/*
	对 circles 中的 每个 圆
	将其绘制在 dstImg 中

	for circle in cirlces:
		invoke cv::cirlce function to draw circle on dstImg
	*/
	void draw_circle(Mat& dstImg, const Circle& circle)
	{
		cv::circle(dstImg, circle.center(), circle.radius, Scalar(0, 100, 255), 2);
	}

	void draw_circle(Mat& dstImg, const vector<Circle>& circles)
	{
		for (auto& circle : circles)
			cv::circle(dstImg, circle.center(), circle.radius, Scalar(0, 100, 255), 2);
	}

	void draw_circle(Mat& dstImg, vector<cv::Vec3f>& circles)
	{
		if (circles.size()) {
			sort(circles.begin(), circles.end(), [](cv::Vec3f& a, cv::Vec3f& b) ->bool {return a[2] < b[2]; });
			size_t i = 0;
			int centerX = cvRound(circles.at(i)[0]);
			int centerY = cvRound(circles.at(i)[1]);
			int radius = cvRound(circles.at(i)[2]);
			MyPoint circleCenter = MyPoint(centerX, centerY);
			cv::circle(dstImg, circleCenter, radius, Scalar(0, 100, 255), 2);
		}
	}

	/*
	通过 borders中给定的元的边界点 拟合出圆
	use a set of points of circle to calculate/simulate the cirlce
	*/
	Circle circle_least_fit(const vector<MyPoint>& border)
	{
		//最小乘二法拟合圆
		double center_x = 0.0f;
		double center_y = 0.0f;
		double radius = 0.0f;

		double sum_x = 0.0f, sum_y = 0.0f;
		double sum_x2 = 0.0f, sum_y2 = 0.0f;
		double sum_x3 = 0.0f, sum_y3 = 0.0f;
		double sum_xy = 0.0f, sum_x1y2 = 0.0f, sum_x2y1 = 0.0f;

		int N = border.size();
		for (int i = 0; i < N; i++)
		{
			double x = border[i].x;
			double y = border[i].y;
			double x2 = x * x;
			double y2 = y * y;
			sum_x += x;
			sum_y += y;
			sum_x2 += x2;
			sum_y2 += y2;
			sum_x3 += x2 * x;
			sum_y3 += y2 * y;
			sum_xy += x * y;
			sum_x1y2 += x * y2;
			sum_x2y1 += x2 * y;
		}

		double C, D, E, G, H;
		double a, b, c;

		C = N * sum_x2 - sum_x * sum_x;
		D = N * sum_xy - sum_x * sum_y;
		E = N * sum_x3 + N * sum_x1y2 - (sum_x2 + sum_y2) * sum_x;
		G = N * sum_y2 - sum_y * sum_y;
		H = N * sum_x2y1 + N * sum_y3 - (sum_x2 + sum_y2) * sum_y;
		a = (H * D - E * G) / (C * G - D * D);
		b = (H * C - E * D) / (D * D - G * C);
		c = -(a * sum_x + b * sum_y + sum_x2 + sum_y2) / N;

		center_x = a / (-2);
		center_y = b / (-2);
		radius = sqrt(a * a + b * b - 4 * c) / 2;

		Circle circle;

		circle.centerX = center_x;
		circle.centerY = center_y;
		circle.radius = radius;
		return move(circle);
	}

	/*
	给一张图片 找到 瞳孔内一点
	寻找策略:
	(1) 消除图像边界： 使图像边界值为 255
	(2) 计算 每9列的 像素和， 找出最小的9 * img.rows
	(3) 在最小的 9 * img.rows 中找出 累加和最小的 9 * 9 区域
	(4) 返回这个 9 * 9 区域的中心点
	*/
	MyPoint find_pupil_inner_point(const Mat& img)
	{
		Mat dealingImg = img.clone();

		eliminate_border(dealingImg, dealingImg.cols / 6, dealingImg.rows / 6);

		int sumCol = 9;
		int sumRow = dealingImg.rows;
		int minSumCol = sumCol * sumRow * 255;
		int minSumRow = minSumCol;
		int minColStart;
		int minRowStart;

		int sumVal;
		for (int i = 0; i < dealingImg.cols - sumCol; i += 1)
		{
			sumVal = accumulate(dealingImg, i, 0, sumCol, sumRow);
			if (sumVal < minSumCol)
			{
				minSumCol = sumVal;
				minColStart = i;
			}
		}
		for (int j = 0; j < dealingImg.rows - sumCol; ++j)
		{
			sumVal = accumulate(dealingImg, minColStart, j, 9, 9);
			if (sumVal < minSumRow)
			{
				minSumRow = sumVal;
				minRowStart = j;
			}
		}
		return MyPoint(minColStart + sumCol / 2, minRowStart + sumCol / 2);
	}


	/*
	将图像的像素点值 按 size分布，没size个像素点的值 为这个区域的像素点平均值
	use average pixel value of every area=size in src to fill each area in dst
	*/
	void average_pixel(Mat& src, Mat& dst, Size size = Size(4, 4))
	{
		for (int i = 0; i < src.cols - size.width; i += size.width)
		{
			for (int j = 0; j < src.rows - size.height; j += size.height)
			{
				int sum = accumulate(src, i, j, size.width, size.height);
				int average = sum / (size.width * size.height);
				fill_value(dst, average, i, j, size.width, size.height);
			}
		}
	}

	/*
	种子生长,
	seed 为瞳孔内一点，maxval 为 判定为瞳孔内一点的像素最大值
	if (img.at(point) < maxval) 则 point 在瞳孔内
	*/
	Mat seed_grow(const Mat& img, MyPoint seed, uchar maxval = 100) {

		std::set<MyPoint> xy_temp;
		std::set<MyPoint> xy;

		xy_temp.insert(seed);

		int xyadd[8][2] = {
			{ -1,-1 },
			{ -1,0 },
			{ -1,1 },
			{ 0,-1 },
			{ 0,1 },
			{ 1,-1 },
			{ 1,0 },
			{ 1,1 }
		};


		while (!xy_temp.empty())
		{
			MyPoint temp = *(xy_temp.begin());
			for (int i = 0; i < 8; i++)
			{
				int col = temp.x + xyadd[i][1];
				int row = temp.y + xyadd[i][0];
				col = col > 0 ? col : 0;
				row = row > 0 ? row : 0;

				if (img.at<uchar>(row, col) < maxval)
				{
					MyPoint&& p = MyPoint(col, row);
					if (xy.count(p) == 0)
					{
						xy_temp.insert(p);
					}
				}
			}
			xy.insert(temp);
			xy_temp.erase(temp);
		}

		cv::Mat bitImg(img.rows, img.cols, CV_8UC1, cv::Scalar(0));

		for (auto& c : xy)
		{
			bitImg.at<uchar>(c.y, c.x) = 255;
		}

		return bitImg;
	}

	/*
	通过给定一个 2值图(0 或 255) 来找出瞳孔的边界点
	*/
	vector<MyPoint> get_border_points(Mat& bitImg)
	{
		vector<MyPoint> border;
		int row = bitImg.rows;
		int col = bitImg.cols;

		int arr_in[] = {
			row * 4 / 16,
			row * 5 / 16,
			row * 6 / 16,
			row * 7 / 16,
			row * 8 / 16,
			row * 9 / 16,
			row * 10 / 16,
			row * 11 / 16,
			row * 12 / 16
		};
		int arr_jn[] = {
			col * 4 / 16,
			col * 5 / 16,
			col * 6 / 16,
			col * 7 / 16,
			col * 8 / 16,
			col * 9 / 16,
			col * 10 / 16,
			col * 11 / 16,
			col * 12 / 16
		};


		vector<int> in(arr_in, arr_in + SIZE(arr_in));
		vector<int> jn(arr_jn, arr_jn + SIZE(arr_jn));


		int j0, jt;

		for (int ini = 0; ini < in.size(); ini++)
		{
			j0 = 0;
			jt = -1;
			for (int j = 0; j < bitImg.cols; j++)
			{
				if (bitImg.at<uchar>(in[ini], j) == 255)
				{
					if (jt == -1)
					{
						j0 = j;
						jt = j;
					}
					else
					{
						jt = j;
					}
				}
			}

			if (j0 != 0)
			{
				border.push_back(MyPoint(j0, in[ini]));
				border.push_back(MyPoint(jt, in[ini]));

			}
		}

		for (int jni = 0; jni < jn.size(); jni++)
		{
			j0 = 0;
			jt = -1;
			for (int i = 0; i < bitImg.rows; i++)
			{
				if (bitImg.at<uchar>(i, jn[jni]) == 255)
				{
					if (jt == -1)
					{
						j0 = i;
						jt = i;
					}
					else
					{
						jt = i;
					}
				}
			}

			if (j0 != 0)
			{
				border.push_back(MyPoint(jn[jni], j0));
				border.push_back(MyPoint(jn[jni], jt));
			}
		}
		return move(border);
	}

}

// 简写 customize 为 cs
namespace cs = customize;


