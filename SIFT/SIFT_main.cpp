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
	cout << "��O�ɶ��� :" << (end - start)/1000 <<  "��" << endl;
	//getchar();//�Ȱ�
	return 0;
}