/*
//bool Iris::test_way_chunk()
//{
//	try
//	{
//		Mat dealingImg = grayImage.clone();
//		int left = -1;
//		int right = -1;
//		int up = -1;
//		int down = -1;
//
//		eliminate_border(dealingImg, dealingImg.cols / 6, dealingImg.rows / 6);
//
//		MyPoint firstPoint = cs::find_pupil_inner_point(dealingImg);
//
//		cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE,
//			cv::Size(2 * g_nStructElementsize + 1, 2 * g_nStructElementsize + 1),
//			cv::Point(g_nStructElementsize, g_nStructElementsize));
//		morphologyEx(dealingImg, dealingImg, cv::MORPH_CLOSE, element);
//		morphologyEx(dealingImg, dealingImg, cv::MORPH_OPEN, element);
//
//		//unsigned char beforehand = getBeforehand() + 30;
//		int beforehand = 500000;
//
//		for (int i = firstPoint.x; i > 0; --i)
//		{
//			if (rectThresholdCount(dealingImg, MyPoint(i, firstPoint.y), 10) > beforehand)
//			{
//				left = i;
//				break;
//			}
//		}
//		for (int i = firstPoint.x; i < dealingImg.cols; ++i)
//		{
//			if (rectThresholdCount(dealingImg, MyPoint(i, firstPoint.y), 10) > beforehand)
//			{
//				right = i;
//				break;
//			}
//		}
//		firstPoint.x = left + ((right - left) >> 1);
//		
//		for (int i = firstPoint.y; i > 0; --i)
//		{
//			if (rectThresholdCount(dealingImg, MyPoint(firstPoint.x, i), 10) > beforehand)
//			{
//				up = i;
//				break;
//			}
//		}
//		for (int i = firstPoint.y; i > dealingImg.rows; ++i)
//		{
//			if (rectThresholdCount(dealingImg, MyPoint(firstPoint.x, i), 10) > beforehand)
//			{
//				down = i;
//				break;
//			}
//		}
//		firstPoint.y = up + ((down - up) >> 1);
//		pupilCircle.centerX = firstPoint.x;
//		pupilCircle.centerY = firstPoint.y;
//		pupilCircle.radius = 10;
//	}
//	catch (const exception& e)
//	{
//		cout << getImageName() << endl;
//		cerr << e.what() << endl;
//		return false;
//	}
//	return true;
//}
//
//int Iris::rectThresholdCount(Mat img, MyPoint centre, int size)
//{
//	int count = 0;
//	int _size = size * 2 + 1;
//	for (int i = centre.x - size; i <= centre.x + size; ++i)
//	{
//		uchar *pData = img.ptr(i);
//		for (int j = centre.y - size; j <= centre.y + size; ++j)
//		{
//			count += (int)pData[j];
//		}
//	}
//
//	return count;
//}

//unsigned char Iris::rectThresholdCount(MyPoint centre, int size)
//{
//	int count = 0;
//	int _size = size * 2 + 1;
//	for (int i = centre.x - size; i <= centre.x + size; ++i)
//	{
//		uchar *pData = grayImage.ptr(i);
//		for (int j = centre.y - size; j <= centre.y + size; ++j)
//		{
//			count += (int)pData[j];
//		}
//	}
//	
//	return (unsigned char)(count / (_size * _size));
//}

//unsigned char Iris::getBeforehand()
//{
//	unsigned char nThreshold = 0;
//	// 计算直方图 histogram
//	int fHistogram[256] = { 0 };
//	for (int y = 0; y < grayImage.rows; ++y)
//	{
//		uchar *pData = grayImage.ptr(y);
//		for (int x = 0; x < grayImage.cols; ++x)
//		{
//			fHistogram[*pData++]++;
//		}
//	}
//	int max = 0;
//	for (int i = 0; i < 100; ++i)
//	{
//		if (fHistogram[i] > max)
//		{
//			max = fHistogram[i];
//			nThreshold = (unsigned char) i;
//		}
//	}
//	return nThreshold;
//}
*/

