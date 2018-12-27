#pragma once
// 定义CASIA-Iris-Interval的文件路径
#define BASE_DIR "F://CASIA-Iris-Interval"

// 重定向标准错误流文件
#define ERROR_FILE "C://Users//Otto//Desktop//error_info.txt"

// 定义从哪个文件开始读和输出图像 (CASIA-Iris-Interval 内部的文件从 001 到 249, 所以RANGE_LEFT 和 RANGE_RIGHT 的值 可以从 1到249
// 建议 RANGE_RIGHT - RANGE_LEFT <=5， 然后每次运行完之后再增加五 
/*
以 要测 前30组为例
第一次:
#define RANGE_LEFT (1)
#define RANGE_RIGHT (5)
.... 运行, 记录结果

第二次:
#define RANGE_LEFT (6)
#define RANGE_RIGHT (10)
.... 运行, 记录结果
以此类推
*/

// 线程数
#define CORES 8

#define RANGE_LEFT (231)

#define RANGE_RIGHT (249)

// 外圆与内圆同心 / 不同心 
#define DIFF_CENTER // 不同心
// #define SAME_CENTER // 同心

// 积分用几圈 : default = 1
#define VERBOSE // 3圈