#include"SIFT.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<Windows.h>
#define PI 3.14159265358979323846
#define E  2.718281828459045235
/* test */
using namespace std;

struct RGBDATA {
	BYTE rgbtRed;
	BYTE rgbtGreen;
	BYTE rgbtBlue;
};

SIFT::SIFT(char* inname)
{
	Readbmp(inname);
	RGB2Gray();
}

SIFT::~SIFT()
{
	for (int i = 0; i < Height; ++i)
	{
		delete[] color[i];
	}
	delete[] color;
}

void SIFT::Readbmp(char * inname)
{
	int fix;
	/* 開啟檔案 */
	fstream fp;
	fp.open(inname, ios::in | ios::binary);
	/* 先讀取檔頭資訊 */
	fp.read((char*)&FileHeader, sizeof(BITMAPFILEHEADER));
	fp.read((char*)&InfoHeader, sizeof(BITMAPINFOHEADER));

	/* 取得圖寬及圖高 */
	Width = InfoHeader.biWidth;
	Height = InfoHeader.biHeight;
	
	color = new RGBTRIPLE*[Height];
	for (int i = 0; i < Height; ++i)
	{
		color[i] = new RGBTRIPLE[Width];
	}

	cout << "Width : " << Width << endl;
	cout << "Height : " << Height << endl;

	/* 計算每列需略過的 bytes 數 */
	if ((Width * 3) % 4 != 0)
		fix = 4 - ((Width * 3) % 4);
	else
		fix = 0;

	/* 圖入各像素資訊 */
	for (int j = Height - 1; j >= 0; --j)
	{
		for (int i = 0; i<Width; ++i)
		{
			RGBDATA rgb;
			fp.read((char*)&rgb, sizeof(RGBTRIPLE));
			color[j][i].rgbtBlue = rgb.rgbtBlue;
			color[j][i].rgbtGreen = rgb.rgbtGreen;
			color[j][i].rgbtRed = rgb.rgbtRed;
		}
		/* 略過各列多餘的資訊 */
		BYTE   ByteBuf;
		for (int n = 0; n<fix; ++n) {
			fp.read((char*)&ByteBuf, sizeof(BYTE));
		}
	}

	//****************************輸出
	OutBMP();
	OutRAW();
	//**************************
}
//*** 檢查是否為可用的特徵點 ***//
void SIFT::FeatureCheck()
{
	
}
//*** 擷取特徵點 ***//
void SIFT::getFeature(float** kaidanImage, int InWidth, int InHeight)
{
	for (int k = 1; k < 2; k++)
	{
		for (int j = 1; j < (InHeight - 1); j++)
		{
			for (int i = 1; i < (InWidth - 1); i++)
			{
				bool rof = checkmaxmin(kaidanImage, i, j, k,InWidth);
				if (rof)
				{
					//處理此特徵點
					float xxx = kaidanImage[1][j*InWidth+i];
				}
			}
		}
	}
}
bool SIFT::checkmaxmin(float ** doImage, int nowx, int nowy, int nowz,int InWidth)//未完待續
{
	float value = doImage[nowz][nowy*InWidth + nowx];
	bool maxc, minc;
	maxc = minc = true;
	for (int k = (nowz - 1); k <= (nowz + 1); k++)
	{
		for (int j = (nowy - 1); j <= (nowy + 1); j++)
		{
			for (int i = (nowx - 1); i < (nowx + 1); i++)
			{
				if (doImage[k][j*InWidth + i] > value)
					maxc = false;
				if (doImage[k][j*InWidth + i] < value)
					minc = false;
			}
		}
	}
	if (!(maxc | minc))
		return false;
	else
		return true;
}
//*** 高斯模糊 ***//
void SIFT::GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma)
{
	float a, b, c,total;
	a = pow(E, -2 / (2 * sigma*sigma)) / (2 * PI*sigma*sigma);
	b = pow(E, -1 / (2 * sigma*sigma)) / (2 * PI*sigma*sigma);
	c = 1 / (2 * PI*sigma*sigma);
	total = 4 * a + 4 * b + c;
	a /= total;
	b /= total;
	c /= total;

	float mask[9] = { a,b,a,b,c,b,a,b,a };
	float block[9];
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (i == 0 && j == 0)
			{
				block[0] = 0;
				block[1] = 0;
				block[2] = 0;

				block[3] = 0;
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = doImage[inz][(j + 0)*InWidth + (i + 1)];

				block[6] = 0;
				block[7] = doImage[inz][(j + 1)*InWidth + (i + 0)];
				block[8] = doImage[inz][(j + 1)*InWidth + (i + 1)];
			}
			else if (i == 0 && j == (InHeight - 1))
			{
				block[0] = 0;
				block[1] = doImage[inz][(j - 1)*InWidth + (i + 0)];
				block[2] = doImage[inz][(j - 1)*InWidth + (i + 1)];

				block[3] = 0;
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = doImage[inz][(j + 0)*InWidth + (i + 1)];

				block[6] = 0;
				block[7] = 0;
				block[8] = 0;
			}
			else if (i == (InWidth - 1) && j == 0)
			{
				block[0] = 0;
				block[1] = 0;
				block[2] = 0;

				block[3] = doImage[inz][(j + 0)*InWidth + (i - 1)];
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = 0;

				block[6] = doImage[inz][(j + 1)*InWidth + (i - 1)];
				block[7] = doImage[inz][(j + 1)*InWidth + (i + 0)];
				block[8] = 0;
			}
			else if (i == (InWidth - 1) && j == (InHeight - 1))
			{
				block[0] = doImage[inz][(j - 1)*InWidth + (i - 1)];
				block[1] = doImage[inz][(j - 1)*InWidth + (i + 0)];
				block[2] = 0;

				block[3] = doImage[inz][(j + 0)*InWidth + (i - 1)];
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = 0;

				block[6] = 0;
				block[7] = 0;
				block[8] = 0;
			}
			else if (i == 0)
			{
				block[0] = 0;
				block[1] = doImage[inz][(j - 1)*InWidth + (i + 0)];
				block[2] = doImage[inz][(j - 1)*InWidth + (i + 1)];

				block[3] = 0;
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = doImage[inz][(j + 0)*InWidth + (i + 1)];

				block[6] = 0;
				block[7] = doImage[inz][(j + 1)*InWidth + (i + 0)];
				block[8] = doImage[inz][(j + 1)*InWidth + (i + 1)];
			}
			else if (i == (InWidth - 1))
			{
				block[0] = doImage[inz][(j - 1)*InWidth + (i - 1)];
				block[1] = doImage[inz][(j - 1)*InWidth + (i + 0)];
				block[2] = 0;

				block[3] = doImage[inz][(j + 0)*InWidth + (i - 1)];
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = 0;

				block[6] = doImage[inz][(j + 1)*InWidth + (i - 1)];
				block[7] = doImage[inz][(j + 1)*InWidth + (i + 0)];
				block[8] = 0;
			}
			else if (j == 0)
			{
				block[0] = 0;
				block[1] = 0;
				block[2] = 0;

				block[3] = doImage[inz][(j + 0)*InWidth + (i - 1)];
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = doImage[inz][(j + 0)*InWidth + (i + 1)];

				block[6] = doImage[inz][(j + 1)*InWidth + (i - 1)];
				block[7] = doImage[inz][(j + 1)*InWidth + (i + 0)];
				block[8] = doImage[inz][(j + 1)*InWidth + (i + 1)];
			}
			else if (j == (InHeight - 1))
			{
				block[0] = doImage[inz][(j - 1)*InWidth + (i - 1)];
				block[1] = doImage[inz][(j - 1)*InWidth + (i + 0)];
				block[2] = doImage[inz][(j - 1)*InWidth + (i + 1)];

				block[3] = doImage[inz][(j + 0)*InWidth + (i - 1)];
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = doImage[inz][(j + 0)*InWidth + (i + 1)];

				block[6] = 0;
				block[7] = 0;
				block[8] = 0;
			}
			else
			{
				block[0] = doImage[inz][(j - 1)*InWidth + (i - 1)];
				block[1] = doImage[inz][(j - 1)*InWidth + (i + 0)];
				block[2] = doImage[inz][(j - 1)*InWidth + (i + 1)];

				block[3] = doImage[inz][(j + 0)*InWidth + (i - 1)];
				block[4] = doImage[inz][(j + 0)*InWidth + (i + 0)];
				block[5] = doImage[inz][(j + 0)*InWidth + (i + 1)];

				block[6] = doImage[inz][(j + 1)*InWidth + (i - 1)];
				block[7] = doImage[inz][(j + 1)*InWidth + (i + 0)];
				block[8] = doImage[inz][(j + 1)*InWidth + (i + 1)];
			}

			for (int v = 0; v < 9; ++v)
			{
				float add = 0;
				 add += block[v] * mask[v];
				 doImage[inz][j*InWidth + i] = add;
			}
		}
	}
}
//*** 高斯差 ***//
void SIFT::GusC(float** FinalImage,int inz, int InWidth, int InHeight)
{
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			float value1 = FinalImage[inz][j*InWidth + i];
			float value2 = FinalImage[inz+1][j*InWidth + i];
			value1 -= value2;
			FinalImage[inz][j*InWidth + i] = value1;
		}
	}
}

void SIFT::Feature()
{
	//************** 以圖片大小計算出能把金字塔做多高 ****************//
	float SizeMin;
	int n = 1;
	float Limit = 1;
	while (Limit < Width || Limit < Height)
	{
		Limit = pow(2.0, n);
		++n;
	}
	if (n < 7) { n = 1; }
	else { n -= 6; }
	Limit = pow(0.5, n);
	//********************* 放大縮小以及高斯模糊 *************************
	for (int size = 1; size >= Limit; size/=2)
	{
		float** KaIDaN;
		KaIDaN = new float*[4];
		int inWidth = Width*size;
		int inHeight = Height*size;
		//高斯模糊
		for (int i = 0; i < 5; ++i)
		{
			KaIDaN[i] = new float[inWidth*inHeight];
			ZoomInOut(KaIDaN[i], inWidth, inHeight);
			GusB(KaIDaN, i, inWidth, inHeight);
		}
		//高斯差
		for (int i = 0; i < 4; i++)
		{
			GusC(KaIDaN, i, inWidth, inHeight);
		}
		//擷取特徵點
		getFeature(KaIDaN, inWidth, inHeight);
		for (int i = 0; i < 5; i++)
		{
			delete[] KaIDaN[i];
		}
	}
}
//*** 放大縮小 ***//
void SIFT::ZoomInOut(float* doImage, int InWidth, int InHeight)
{
	float XX1, XX2, XXYY;
	float dx, dy;
	float Io, Jo;
	for (int j = 0; j < InHeight-1; ++j)
	{
		Jo = j * (Height - 1) / (InHeight - 1);
		dy = 1 - ((j * (Height - 1) / (InHeight - 1)) - Jo);
		for (int i = 0; i < InWidth-1; ++i)
		{
			Io = i * (Width - 1) / (InWidth - 1);
			dx = 1 - ((i * (Width - 1) / (InWidth - 1)) - Io);
			XX1 = gray[(Io + 0) + (Jo + 0)*Width] * dx + gray[(Io + 1) + (Jo + 0)*Width] * (1 - dx);
			XX2 = gray[(Io + 0) + (Jo + 1)*Width] * dx + gray[(Io + 1) + (Jo + 1)*Width] * (1 - dx);
			XXYY = XX1*dy + XX2*(1 - dy);
			doImage[j*InWidth + i] = XXYY;
		}
		XX1 = gray[(Jo + 1)*Width - 1];
		XX2 = gray[(Jo + 2)*Width - 1];
		XXYY = XX1*dy + XX2*(1 - dy);
		doImage[(j + 1)*InWidth - 1] = XXYY;
	}
	for (int i = 0; i < InWidth-1; i++)
	{
		Io = i * (Width - 1) / (InWidth - 1);
		dx = 1 - ((i * (Width - 1) / (InWidth - 1)) - Io);
		XX1 = gray[(Height - 1)*Width + (Io + 0)];
		XX2 = gray[(Height - 1)*Width + (Io + 1)];
		XXYY = XX1*dx + XX2*(1 - dx);
		doImage[(InHeight - 1)*InWidth + i] = XXYY;
	}
	doImage[InHeight*InWidth - 1] = gray[Height*Width - 1];
}
//*** 轉灰階 ***//
void SIFT::RGB2Gray()
{
	gray.resize(Height*Width);
	for (int j = 0; j < Height; ++j)
	{
		for (int i = 0; i < Width; ++i)
		{
			gray[j*Width + i] = color[j][i].rgbtRed*0.299 + color[j][i].rgbtGreen*0.587 + color[j][i].rgbtBlue*0.114;
			gray[j*Width + i] /= 255;
		}
	}
}
//*** 將檔案輸出成BMP格式 ***//
void SIFT::OutBMP()
{
	int fix;
	/* 計算每列需略過的 bytes 數 */
	if ((Width * 3) % 4 != 0)
		fix = 4 - ((Width * 3) % 4);
	else
		fix = 0;

	fstream out;
	out.open("out.bmp", ios::out | ios::binary);
	out.write((char*)&FileHeader, sizeof(BITMAPFILEHEADER));
	out.write((char*)&InfoHeader, sizeof(BITMAPINFOHEADER));

	for (int j = Height-1; j >= 0; --j)
	{
		for (int i = 0; i < Width; ++i)
		{
			RGBDATA rgb;
			rgb.rgbtBlue = color[j][i].rgbtBlue;
			rgb.rgbtGreen = color[j][i].rgbtGreen;
			rgb.rgbtRed = color[j][i].rgbtRed;
			out.write((char*)&rgb, sizeof(RGBTRIPLE));
		}
		/* 略過各列多餘的資訊 */
		for (int n = 0; n<fix; ++n) {
			BYTE   ByteBuf;
			out.write((char*)&ByteBuf, sizeof(ByteBuf));
		}
	}
}
//*** 將檔案輸出成RAW格式 ***//
void SIFT::OutRAW()
{
	fstream out;
	out.open("out.raw", ios::out | ios::binary);

	for (int j = 0; j < Height; ++j)
	{
		for (int i = 0; i < Width; ++i)
		{
			out.write((char*)&color[j][i], sizeof(RGBTRIPLE));
		}
	}
}
