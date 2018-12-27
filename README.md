# 代码说明
该代码的功能实现了对虹膜的分割操作，使用的训练样本图片包含在CASIA-Iris-Interval中，请在使用时更改interface.h里面的文件路径（BASE_DIR)到你下载CASIA-Iris-Interval的
路径上。

# 运行要求
该代码因为 include 了 Windows.h 故只能在Windows环境下运行，同时代码中使用了C++11,14,17的语言标准，所以请在运行时在工程属性里设置语言标准为C++17或C
++最新标准，同时关掉SDL检查。另外的，因为代码大量的调用了opencv的库函数，所以运行时必须安装好opencv环境，配置好环境变量后方可运行。

# CASIA-Iris-Interval
中科院的公开的虹膜图像库，请将CASIA-Iris-Interval.0.tar.gz和CASIA-Iris-Interval.1.gz和并解压缩到CASIA-Iris-Interval的文件夹
