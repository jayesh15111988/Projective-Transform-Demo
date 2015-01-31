transform : transform.cpp SImage.h SImageIO.h DTwoDimArray.h
	g++ -O3 -o transform transform.cpp -I . -lpng
