#pragma once

#include "irishandler.h"

class Iris : public IrisHandler {
	/**
	public 继承 IrisHandler 实现 getPupilCirlce 方法
	*/
protected:
	Mat originalImage;

public:
	Iris() = default;

	Iris(std::string imagePath);

	bool getPupilCircle() noexcept override;

	bool getIrisCircle() noexcept override;

	void add_image(std::string imagePath);

	bool original_find_pupilCircle();

	bool find_pupilCircle();

	// 圆心区域左上角点 和 右下角点

	pair<MyPoint, MyPoint> getPupilCenterArea(int xRadius = 8, int yRadius = 8);

	bool findPupilEdge();

	MyPoint findPupilPoint();

	int calCircleSum(cv::Mat &img, int x, int y, int r);

	int calCircleSum_pupil(cv::Mat& img, int x, int y, int r);

	bool findIrisEdge(const cv::Mat &src);

	void drawCircle() noexcept;
};
