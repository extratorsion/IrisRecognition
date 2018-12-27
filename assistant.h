#pragma once

#include "base.h"
/*
继承自cv::Point，与cv::Point可以混用
*/
class MyPoint : public cv::Point_<int>
{
public:

	using Point_<int>::Point_;

	MyPoint() = default;

	MyPoint(const MyPoint& op) = default;

	MyPoint(const cv::Point& op)
	{
		x = op.x, y = op.y;
	}

	operator cv::Point() const 
	{
		return cv::Point(x, y); 
	}

	void set(int x, int y) 
	{
		this->x = x, this->y = y;
	}

	const MyPoint& operator=(const MyPoint& op)
	{
		this->x = op.x, this->y = op.y;
		return *this;
	}

	const MyPoint& operator=(const Point& op)
	{
		this->x = op.x, this->y = op.y;
		return *this;
	}

	bool operator<(const Point& op) const
	{
		return x < op.x ? true : x > op.x ? false : y >= op.y ? false : true;
	}

	bool operator>(const Point& op) const
	{
		return x < op.x ? false : x > op.x ? true : y >= op.y ? true : false;
	}

	bool operator==(const Point& op) const
	{
		return x == op.x && y == op.y;
	}

	bool operator!=(const Point& op) const
	{
		return !operator==(op);
	}

	// 计算该点相聚另一点的距离，返回int
	int apart(const Point& op) const
	{
		return cvRound(sqrt(
			pow(op.x - x, 2) + pow(op.y - y, 2)
		));
	}

	friend ostream& operator<<(ostream& os, MyPoint& point)
	{
		os << "MyPoint<" << point.x << ", " << point.y << ">"; return os;
	}

	// 返回 MyPoint(this->x + xdraft, this->y);
	MyPoint Xdraft(int xdraft)
	{
		return MyPoint(this->x + xdraft, this->y);
	}

	// 返回 MyPoint(this->x, this->y + ydraft)
	MyPoint Ydraft(int ydraft)
	{
		return MyPoint(this->x, this->y + ydraft);
	}

	// 返回 MyPoint(this->x +  draft, this->y + draft)
	MyPoint draft(int draft) 
	{
		return MyPoint(this->x +  draft, this->y + draft);
	}

	// 返回 MyPoint(this->x +  xdraft, this->y + ydraft)
	MyPoint draft(int xdraft, int ydraft)
	{
		return MyPoint(this->x +  xdraft, this->y + ydraft);
	}
};

template <>
struct hash<MyPoint> {
	size_t operator()(const MyPoint& p) {
		return size_t(p.x * p.y + (p.x ^ p.y));
	}
};


class Circle {
public:
	int centerX = -1;
	int centerY = -1;
	int radius = -1;

public:

	Circle() = default;

	Circle(const Circle& circle) = default;

	Circle(int x, int y, int radius)
		: centerX(x), centerY(y), radius(radius)
	{}

	Circle(const cv::Vec3i& veci)
		: centerX(veci[0]), centerY(veci[1]), radius(veci[2])
	{}

	Circle(const cv::Vec3f& vecf)
		: centerX(cvRound(vecf[0])), centerY(cvRound(vecf[1])), radius(cvRound(vecf[2]))
	{}

	Circle(const cv::Vec3d& vecd)
		: centerX(cvRound(vecd[0])), centerY(cvRound(vecd[1])), radius(cvRound(vecd[2]))
	{}

	MyPoint center() const { return MyPoint(centerX, centerY); }

	bool is_assigned() const { return radius != -1 && centerX != -1 && centerY != -1; }

	bool is_valid() const { return radius > 0 && centerX >= 0 && centerY >= 0; }

	bool in_cirlce(const Point& point) const {
		double distance = sqrt((pow((point.x - centerX), 2) + pow((point.y - centerY), 2)));
		return distance < ((double)radius);
	}

	bool in_circle(int row, int col) const {
		int distance = cvRound(sqrt(pow(col - centerX, 2) + pow(row - centerY, 2)));
		return distance < radius;
	}

	const Circle& operator=(const Circle& circle)
	{
		this->centerX = circle.centerX;
		this->centerY = circle.centerY;
		this->radius = circle.radius;
		return *this;
	}

	friend ostream& operator << (ostream& os, Circle& circle)
	{
		if (circle.is_assigned())
			os << "Circle<" << circle.centerX << ", " << circle.centerY << " |" << circle.radius << ">";
		else
			os << "Circle<" << "UnSigned" << ">";
		return os;
	}

public:
	class CircleWrongArgument : public exception {
		using exception::exception;
		virtual const char* what() const throw() {
			return "CircleWrongArgument: the circle's argument is wrong";
		}
	};
};

template<>
struct hash<Circle> {
	size_t operator()(const Circle& c) {
		return size_t(c.centerX * c.centerY + c.radius);
	}
};