#include<iostream>
#include<ctime>
#include<string>
#include<vector>
#include"SIFT.h"
#define UsePicture "ball_01.bmp"

using namespace std;



int main(int argc, char** argv)
{
	clock_t start,end;
	start = clock();
	//*****************************************************************************
	SIFT ImageSIFT(UsePicture);
	ImageSIFT.Feature();
	/* test */
	//*****************************************************************************
	end = clock();
	cout << "花費時間為 :" << (end - start)/1000 <<  "秒" << endl;
	//getchar();//暫停
	return 0;
}