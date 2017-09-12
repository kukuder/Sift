#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include "SIFT.h"
#define UsePicture1 "kanna.bmp"
#define UsePicture2 "ball_02.bmp"
#define OutPicture1 "Featurename1"
#define OutPicture2 "Featurename2"

using namespace std;



int main(int argc, char** argv)
{
	clock_t start,end;
	start = clock();
	//*****************************************************************************
	SIFT ImageSIFT1(UsePicture1);
	ImageSIFT1.display();
	ImageSIFT1.OutBMP(OutPicture1);

	//SIFT ImageSIFT2(UsePicture2);
	//ImageSIFT2.display();
	//ImageSIFT2.OutBMP(OutPicture2);

	//*****比較特徵點
	//**********************************
	// 1. 新增圖片物件，創造的新圖為原圖等寬、兩倍長
	//    (4.)讀取兩張圖片資料
	// 2. 試著將新圖寫出BMP檔(注意4位元差補問題)
	// 3. 比較特徵描述子(研究資料結構傳輸方面的問題)
	// 4. 將匹配的點相連接
	//**********************************
	//*****************************************************************************
	end = clock();
	cout << "花費時間為 :" << (end - start)/1000 <<  "秒" << endl;
	//getchar();//暫停
	return 0;
}