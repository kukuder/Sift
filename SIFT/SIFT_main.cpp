#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include "SIFT.h"
#define UsePicture1 "dog.bmp"
#define UsePicture2 "dog.bmp"
#define OutPicture1 "Featurename1"
#define OutPicture2 "Featurename2"

using namespace std;



int main(int argc, char** argv)
{
	clock_t start,end;
	start = clock();
	//*****************************************************************************
	SIFT ImageSIFT1(UsePicture1);
	SIFT ImageSIFT2(UsePicture2);
	
	Stitching bigImageSIFT(ImageSIFT1.GetFileHeader(), ImageSIFT1.GetInfoHeader(), ImageSIFT1.Getcolor(), ImageSIFT2.Getcolor(), ImageSIFT1.GetFeatureptr(), ImageSIFT2.GetFeatureptr());//取得標頭檔以及圖片內容
	bigImageSIFT.Check();
	bigImageSIFT.OutBMP("bigpicture");

	ImageSIFT1.display();
	ImageSIFT2.display();
	ImageSIFT1.OutBMP(OutPicture1);
	ImageSIFT2.OutBMP(OutPicture2);
	
	//*****比較特徵點
	//**********************************
	// 4. 將匹配的點相連接
	//**********************************
	//*****************************************************************************
	end = clock();
	cout << "花費時間為 :" << (end - start)/1000 <<  "秒" << endl;
	//getchar();//暫停
	return 0;
}