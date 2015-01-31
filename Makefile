stereo : stereo.cpp SImage.h SImageIO.h DTwoDimArray.h
	g++ -O3 -o stereo stereo.cpp -I . -lpng
