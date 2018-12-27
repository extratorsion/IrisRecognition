#pragma once

#include "assistant.h"

class IrisHandler {

protected:
	Circle pupilCircle;	// 虹膜内圆
	Circle irisCircle;	// 虹膜外圆
	Mat grayImage;        // 单通道灰度图片
	string imagePath;	// 图片的绝对路径

public:

	virtual bool getPupilCircle() noexcept = 0;	// 找到虹膜内圆/瞳孔外圆，将圆的信息存放到pupilCircle 中，失败返回false
	virtual bool getIrisCircle() noexcept = 0; 	//后期的任务 找到虹膜外圆，将圆的信息存放到pupilCircle 中，失败返回false
	string getImageName();
};

