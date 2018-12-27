#pragma once

#include <Windows.h>
#include <thread>
#include "iris.h"
#include "interface.h"
#include <mutex>
#define T

vector<string> filename;

int GetALLgpxFilepathFromfolder(const char * Path)
{
	char szFind[MAX_PATH];
	WIN32_FIND_DATA FindFileData;
	strcpy(szFind, Path);
	strcat(szFind, "//*.*");

	HANDLE hFind = FindFirstFile(szFind, &FindFileData);
	if (INVALID_HANDLE_VALUE == hFind)
		return -1;

	do {

		if (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
		{
			if (strcmp(FindFileData.cFileName, ".") != 0 && strcmp(FindFileData.cFileName, "..") != 0)
			{

				//发现子目录，递归之
				char szFile[MAX_PATH] = { 0 };
				strcpy(szFile, Path);
				strcat(szFile, "//");
				strcat(szFile, FindFileData.cFileName);
				GetALLgpxFilepathFromfolder(szFile);
			}
		}

		else
		{
			//找到文件
			string file = Path;
			file += "//";
			file += FindFileData.cFileName;
			if (*(file.end() - 1) == 'g')
			{
				filename.push_back(file);
			}
		}

	} while (FindNextFile(hFind, &FindFileData));

	FindClose(hFind);

	return 0;
}

inline size_t stringConvertToInteger(const string& str) {
	size_t val = 0;
	size_t len = str.size();
	for (int i = 0; i < len; i++)
		val += pow(10, (len - i - 1)) * (str[i] - '0');
	return val;
}

void run(vector<string> filenames) {
	if (filename.empty())
		return;
	for (auto& _filename: filenames){
		Iris irisObj(_filename);

		if (irisObj.getPupilCircle() == false)
			std::cout << irisObj.getImageName() << " fail in segment pupilCircle" << endl;
#if defined(T)
		else {
			if (irisObj.getIrisCircle() == false)
				std::cout << irisObj.getImageName() << " fail in segment irisCircle" << endl;
			else
				irisObj.drawCircle();
		}
#else 
	irisObj.drawCircle();
#endif
	}
	cv::waitKey();
}

void regular()
{
	cv::namedWindow("Iris Recognition Test");

	int cnt = 0;

	const string path = BASE_DIR;
		
	freopen(ERROR_FILE, "a", stderr); //cerr输出的 文件

	GetALLgpxFilepathFromfolder(path.c_str());		//得到path下所有文件

	const string flag = string("CASIA-Iris-Interval");
	const size_t flagLen = flag.size();

	// 最多线程个数 == CORES 
	vector<thread> threads(CORES);

	vector<vector<string>> filename_list(CORES);

	for (int i = 0; i < filename.size(); i++)
	{
		int ind = filename[i].find(flag);
		string subStr = filename[i].substr(ind + flagLen + 2, 3); // 获得文件名

		int val = stringConvertToInteger(subStr);
		if (val < RANGE_LEFT) {
			continue;
		}
		else if (val >= RANGE_LEFT && val <= RANGE_RIGHT) {
			filename_list[cnt % CORES].push_back(filename[i]);
			++cnt;
		}
		else {
			break;
		}
	}

	// 真正的线程数
	int real_thread_num = cnt > CORES ? CORES : cnt;

	for (int i = 0; i < real_thread_num; ++i) {
		threads.emplace(threads.begin() + i, run, filename_list[i]);
		threads[i].detach(); 
	}

	std::cout << "total " << cnt << " images is testing!" << endl;
}


int main(void)
{
	try {
		regular(); 
	} catch (const exception& e) {
		cout << e.what() << endl;
	}
	system("pause");
	return 0;
}
