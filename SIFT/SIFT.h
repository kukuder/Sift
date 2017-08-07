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
	int x, y;//原圖尺度之未旋轉的座標
	float mm;//強度
	int sita;//包含主方向與負方向的角度
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
	//*** 角點偵測 ***//
	bool Hessian(float* DoImage, int Inx, int Iny, int InWidth);
	//*** 擷取特徵點 ***//
	void getFeature(float** kaidanImage, int InWidth, int InHeight, int insize);
	//*** 檢查該點是否為極值 ***//
	bool checkmaxmin(float** doImage, int nowx, int nowy, int nowz, int InWidth);
	//*** 高斯模糊 ***//
	void GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma = 1.0);
	//*** 計算高斯遮罩的機率(一維) ***//
	void GusPersent(float *inarray, int line,float sigma = 1.0);
	//*** 高斯差 ***//
	void GusC(float** FinalImage, int inz, int InWidth, int InHeight);
	//*** 印出箭頭 ***//
	void display();
	//*** 在特徵點串列後頭加入新的特徵點 ***//
	void AddnewFeaturestruct(int Inx, int Iny, float Inm, int Insita);
	//*** 添加特徵點 ***//
	void AddFeature(float* doImage, int Inx, int Iny, int Inr, int Insize,int InWidth,float sigma);
	//*** 完成一連串尋找特徵點的動作 ***//
	void doFeature();
	//*** 放大縮小 ***//
	void ZoomInOut(float* doImage, int InWidth, int InHeight);
	//*** 轉灰階並且正規化(pixel/255) ***//
	void RGB2Gray();
	//*** 將檔案輸出成BMP格式 ***//
	void OutBMP(std::string outname);
	//*** 將檔案輸出成RAW格式 ***//
	void OutRAW(std::string outname);
};
