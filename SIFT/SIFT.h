#pragma once
#include<Windows.h>
#include<vector>


class SIFT
{
private:
	int Width, Height;
	std::vector <float> gray;
	RGBTRIPLE** color;
	BITMAPFILEHEADER FileHeader;
	BITMAPINFOHEADER InfoHeader;
public:
	SIFT(char* inname);
	~SIFT();
	void Readbmp(char* inname);
	//*** �ˬd�O�_���i�Ϊ��S�x�I ***//
	void FeatureCheck();
	void getFeature(float** kaidanImage, int InWidth, int InHeight);
	//*** �^���S�x�I ***//
	void getFeature(float** upImage);
	//*** �ˬd���I�O�_������ ***//
	bool checkmaxmin(float** doImage, int nowx, int nowy, int nowz, int InWidth);
	//*** �����ҽk ***//
	void GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma = 1.0);
	//*** �����t ***//
	void GusC(float** FinalImage, int inz, int InWidth, int InHeight);
	void Feature();
	//*** ��j�Y�p ***//
	void ZoomInOut(float* doImage, int InWidth, int InHeight);
	//*** ��Ƕ� ***//
	void RGB2Gray();
	//*** �N�ɮ׿�X��BMP�榡 ***//
	void OutBMP();
	//*** �N�ɮ׿�X��RAW�榡 ***//
	void OutRAW();
};
