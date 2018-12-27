
#include "irishandler.h"

string IrisHandler::getImageName() {
	size_t ind = imagePath.rfind("/");
	return imagePath.substr(ind + 1, imagePath.size());
}



