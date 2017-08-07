#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <Windows.h>
#include <string>
/* test */

struct RGBDATA {
	BYTE rgbtRed;
	BYTE rgbtGreen;
	BYTE rgbtBlue;
};

struct Feature
{
	Feature(void) { nextptr = NULL; }
	int x, y;//��Ϥثפ������઺�y��
	float mm;//�j��
	int sita;//�]�t�D��V�P�t��V������
	Feature *nextptr;
};
typedef struct Feature* Featureptr;

class SIFT
{
private:
	
	int Width, Height;
	std::vector <float> gray;
	Featureptr FeatureStart, FeatureEnd;
	RGBTRIPLE** color;
	BITMAPFILEHEADER FileHeader;
	BITMAPINFOHEADER InfoHeader;
public:
	SIFT(char* inname);
	~SIFT();
	void Readbmp(char* inname);
	//*** ���I���� ***//
	bool Hessian(float* DoImage, int Inx, int Iny, int InWidth);
	//*** �^���S�x�I ***//
	void getFeature(float** kaidanImage, int InWidth, int InHeight, int insize);
	//*** �ˬd���I�O�_������ ***//
	bool checkmaxmin(float** doImage, int nowx, int nowy, int nowz, int InWidth);
	//*** �����ҽk ***//
	void GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma = 1.0);
	//*** �p�Ⱚ���B�n�����v(�@��) ***//
	void GusPersent(float *inarray, int line,float sigma = 1.0);
	//*** �����t ***//
	void GusC(float** FinalImage, int inz, int InWidth, int InHeight);
	//*** �L�X�b�Y ***//
	void display();
	//*** �b�S�x�I��C���Y�[�J�s���S�x�I ***//
	void AddnewFeaturestruct(int Inx, int Iny, float Inm, int Insita);
	//*** �K�[�S�x�I ***//
	void AddFeature(float* doImage, int Inx, int Iny, int Inr, int Insize,int InWidth,float sigma);
	//*** �����@�s��M��S�x�I���ʧ@ ***//
	void doFeature();
	//*** ��j�Y�p ***//
	void ZoomInOut(float* doImage, int InWidth, int InHeight);
	//*** ��Ƕ��åB���W��(pixel/255) ***//
	void RGB2Gray();
	//*** �N�ɮ׿�X��BMP�榡 ***//
	void OutBMP(std::string outname);
	//*** �N�ɮ׿�X��RAW�榡 ***//
	void OutRAW(std::string outname);
};
