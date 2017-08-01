#pragma once
#include<Windows.h>
#include<vector>
/* test */

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
	//*** 檢查是否為可用的特徵點 ***//
	void FeatureCheck();
	//*** 角點偵測 ***//
	void Hessian();
	//*** 擷取特徵點 ***//
	void getFeature(float** kaidanImage, int InWidth, int InHeight, int insize);
	//*** 檢查該點是否為極值 ***//
	bool checkmaxmin(float** doImage, int nowx, int nowy, int nowz, int InWidth);
	//*** 高斯模糊 ***//
	void GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma = 1.0);
	//*** 高斯差 ***//
	void GusC(float** FinalImage, int inz, int InWidth, int InHeight);
	void Feature();
	//*** 放大縮小 ***//
	void ZoomInOut(float* doImage, int InWidth, int InHeight);
	//*** 轉灰階 ***//
	void RGB2Gray();
	//*** 將檔案輸出成BMP格式 ***//
	void OutBMP();
	//*** 將檔案輸出成RAW格式 ***//
	void OutRAW();
};
