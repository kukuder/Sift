#pragma once
#include <Windows.h>


struct RGBDATA;

struct Feature;
typedef struct Feature* Featureptr;

class Stitching
{
private:
	int Width, Height;
	RGBTRIPLE** color;
	BITMAPFILEHEADER FileHeader;
	BITMAPINFOHEADER InfoHeader;
public:
	Stitching(BITMAPFILEHEADER Filedata, BITMAPINFOHEADER Infodata);
	~Stitching();
	void OutBMP(std::string outname);
	void OutRAW(std::string outname);
	void WriteImage(RGBTRIPLE** image1, RGBTRIPLE** image2);
};


class SIFT
{
private:
	
	int Width, Height;
	float *gray;
	Featureptr FeatureStart, FeatureEnd, FeatureNow;
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
	void getFeature(float** kaidanImage, float** DogImage, int InWidth, int InHeight, int insize, float sigma);
	//*** 檢查該點是否為極值 ***//
	bool checkmaxmin(float** doImage, int nowx, int nowy, int nowz, int InWidth);
	//*** 高斯模糊 ***//
	void GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma = 1.0);
	//*** 計算高斯遮罩的機率(一維) ***//
	void GusPersent(float *inarray, int line,float sigma = 1.0);
	//*** 高斯差 ***//
	void GusC(float** FinalImage, float** DogImage, int inz, int InWidth, int InHeight);
	//*** 印出箭頭 ***//
	void display();
	//*** 在特徵點串列後頭加入新的特徵點 ***//
	void AddnewFeaturestruct(int Inx, int Iny, int Insize,int kai, int sigmaOCT, float Inm, int Insita);
	//*** 添加特徵點 ***//
	void AddFeature(float* doImage, int Inx, int Iny, int Inr, int Insize,int InWidth,float sigma, int kai);
	//*** 完成一連串尋找特徵點的動作 ***//
	void doFeature();
	//*** 產生特徵描述子 ***//
	void FeatureDescrip(float** kaidaImage);
	//*** 放大縮小 ***//
	void ZoomInOut(float* doImage, int InWidth, int InHeight);
	//*** 轉灰階並且正規化(pixel/255) ***//
	void RGB2Gray();
	//*** 取得正規化的灰階圖片 ***//
	void GetGrayPicture(float* doImage);
	//*** 將檔案輸出成BMP格式 ***//
	void OutBMP(std::string outname);
	//*** 將檔案輸出成RAW格式 ***//
	void OutRAW(std::string outname);
	//*** 鏈結串列解構 ***//
	void Destroy();
	//*** 高斯模糊(通用版) ***//
	void Gusto(float* doImage, int InWidth, int InHeight, float sigma = 1.0);
	//*** 輸出BMP(檢查用) ***//
	void OutBMPCheck(std::string outname, int inwidth, int inheight, float ** doImage);
};


//*** 座標轉換 ***//
void CoordinateChange(int* deltaX, int* deltaY, float sita);
//*** 求m(強度)值 ***//
float getM(float x0, float x1, float y0, float y1);
//*** 求sita值(徑度) ***//
float getSita(float x0, float x1, float y0, float y1);