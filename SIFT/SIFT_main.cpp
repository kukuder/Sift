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

	//*****����S�x�I
	//**********************************
	// 1. �s�W�Ϥ�����A�гy���s�Ϭ���ϵ��e�B�⭿��
	//    (4.)Ū����i�Ϥ����
	// 2. �յ۱N�s�ϼg�XBMP��(�`�N4�줸�t�ɰ��D)
	// 3. ����S�x�y�z�l(��s��Ƶ��c�ǿ�譱�����D)
	// 4. �N�ǰt���I�۳s��
	//**********************************
	//*****************************************************************************
	end = clock();
	cout << "��O�ɶ��� :" << (end - start)/1000 <<  "��" << endl;
	//getchar();//�Ȱ�
	return 0;
}