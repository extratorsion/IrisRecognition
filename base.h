#pragma once

#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\opencv.hpp>
#include <iostream>
#include <functional>
using cv::Mat;
using cv::Mat_;
using cv::Scalar;
using cv::Point;
using cv::Size;
using namespace std;

#ifndef SIZE
#define SIZE(arr) (sizeof(arr)/sizeof(arr[0]))
#endif // !SIZE

#ifndef STR
#define STR(val) #val
#endif // !STR
#ifndef MERGE
#define MGERGE(x, y) x##y
#endif // !MERGE

#ifndef $show
#define $show(expr) std::cout<< STR(expr: )<<expr<<std::endl
#endif 
