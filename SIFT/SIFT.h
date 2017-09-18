#pragma once
#include <Windows.h>


struct RGBDATA;

struct Feature;
typedef struct Feature* Featureptr;

class SIFT
{
private:
	
	int Width, Height;
	double *gray;
	Featureptr FeatureStart, FeatureEnd, FeatureNow;
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
	void getFeature(float** kaidanImage, float** DogImage, int InWidth, int InHeight, float insize, float sigma);
	//*** �ˬd���I�O�_������ ***//
	bool checkmaxmin(float** doImage, int nowx, int nowy, int nowz, int InWidth);
	//*** �����ҽk ***//
	void GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma = 1.0);
	//*** �p�Ⱚ���B�n�����v(�@��) ***//
	void GusPersent(float *inarray, int line,float sigma = 1.0);
	//*** �����t ***//
	void GusC(float** FinalImage, float** DogImage, int inz, int InWidth, int InHeight);
	//*** �L�X�b�Y ***//
	void display();
	//*** �b�S�x�I��C���Y�[�J�s���S�x�I ***//
	void AddnewFeaturestruct(int Inx, int Iny, float Insize,int kai, int sigmaOCT, float Inm, int Insita);
	//*** �K�[�S�x�I ***//
	void AddFeature(float* doImage, int Inx, int Iny, int Inr, float Insize,int InWidth,float sigma, int kai);
	//*** �����@�s��M��S�x�I���ʧ@ ***//
	void doFeature();
	//*** ���ͯS�x�y�z�l ***//
	void FeatureDescrip(float** kaidaImage);
	//*** ��j�Y�p ***//
	void ZoomInOut(float* doImage, int InWidth, int InHeight);
	//*** ��Ƕ��åB���W��(pixel/255) ***//
	void RGB2Gray();
	//*** ���o���W�ƪ��Ƕ��Ϥ� ***//
	void GetGrayPicture(float* doImage);
	//*** �N�ɮ׿�X��BMP�榡 ***//
	void OutBMP(std::string outname);
	//*** �N�ɮ׿�X��RAW�榡 ***//
	void OutRAW(std::string outname);
	//*** �쵲��C�Ѻc ***//
	void Destroy();
	//*** �����ҽk(�q�Ϊ�) ***//
	void Gusto(float* doImage, int InWidth, int InHeight, float sigma = 1.0);
	//*** ��XBMP(�ˬd��) ***//
	void OutBMPCheck(std::string outname, int inwidth, int inheight, float ** doImage);
	//*** ���o���Y��1 ***//
	BITMAPFILEHEADER GetFileHeader(void);
	//*** ���o���Y��2 ***//
	BITMAPINFOHEADER GetInfoHeader(void);
	//*** ���o�Ϥ����e ***//
	RGBTRIPLE** Getcolor(void);
	//*** ���o�S�x�I�Ƶ��c���}�Y ***//
	Featureptr GetFeatureptr(void);
};

class Stitching
{
private:
	int Width, Height;
	Featureptr FeatureStart1, FeatureStart2;
	RGBTRIPLE** color;
	BITMAPFILEHEADER FileHeader;
	BITMAPINFOHEADER InfoHeader;
public:
	//*** �e�ⶵ�ѼƬ����Y�ɡA�ĤT�B�|�����Ϥ���T�A���B�������S�x�I��T ***//
	Stitching(BITMAPFILEHEADER Filedata, BITMAPINFOHEADER Infodata, RGBTRIPLE** image1, RGBTRIPLE** image2, Featureptr inFeatureptr1, Featureptr inFeatureptr2);
	~Stitching();
	void OutBMP(std::string outname);
	void OutRAW(std::string outname);
	void WriteImage(RGBTRIPLE** image1, RGBTRIPLE** image2);
	//*** �ˬd�O�_���ۦP���S�x�y�z�l ***//
	void Check(void);
	//*** �N�a�J�����I�۳s ***//
	void Link(int x1, int y1, int x2, int y2);
	//*** �p���Ӵy�z�l�������ڦ��Z�� ***//
	float EuclideanDistance(std::vector <std::vector <std::vector <float>>> point1, std::vector <std::vector <std::vector <float>>> point2);
};


//*** �y���ഫ ***//
void CoordinateChange(int* deltaX, int* deltaY, float sita);
//*** �Dm(�j��)�� ***//
float getM(float x0, float x1, float y0, float y1);
//*** �Dsita��(�|��) ***//
float getSita(float x0, float x1, float y0, float y1);