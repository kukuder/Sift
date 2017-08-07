#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <Windows.h>
#include <string>
#include "SIFT.h"
#define PI 3.14159265358979323846
#define E  2.718281828459045235
/* test */
using namespace std;



SIFT::SIFT(char* inname)
{
	Readbmp(inname);
	RGB2Gray();
	FeatureStart = new Feature;
	FeatureEnd = FeatureStart;
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
	//OutBMP("out.bmp");
	//OutRAW("out.raw");
	//**************************
}
//*** 角點偵測 ***//
bool SIFT::Hessian(float * DoImage, int Inx, int Iny ,int InWidth)
{
	float Dxx, Dyy, Dxy;
	Dxx = DoImage[(Iny + 0)*InWidth + (Inx + 0)] * 2 - DoImage[(Iny + 0)*InWidth + (Inx - 1)] - DoImage[(Iny + 0)*InWidth + (Inx + 1)];
	Dyy = DoImage[(Iny + 0)*InWidth + (Inx + 0)] * 2 - DoImage[(Iny + 1)*InWidth + (Inx + 0)] - DoImage[(Iny - 1)*InWidth + (Inx + 0)];
	Dxy = DoImage[(Iny + 1)*InWidth + (Inx + 1)] + DoImage[(Iny - 1)*InWidth + (Inx - 1)] - DoImage[(Iny - 1)*InWidth + (Inx + 1)] - DoImage[(Iny + 1)*InWidth + (Inx - 1)];
	Dxy /= 4;

	float Tr, Det;
	Tr = Dxx + Dyy;
	Det = Dxx*Dyy - Dxy*Dxy;
	if ((Tr*Tr / Det) < 12.1)//12.1 是 (r+1)/r  r為10 所得的值
	{
		return true;
	}
	else
	{
		return false;
	}
}
//*** 擷取特徵點 ***//
void SIFT::getFeature(float** kaidanImage, int InWidth, int InHeight, int insize)
{
	for (int k = 1; k < 2; k++)
	{
		int r = insize * 3 * 1.5*(k + 1);
		for (int j = r + 1; j < (InHeight - r - 1); j++)
		{
			for (int i = r + 1; i < (InWidth - r - 1); i++)
			{
				bool rof = checkmaxmin(kaidanImage, i, j, k,InWidth);//極值檢測
				if (rof == true)
				{
					if (Hessian(kaidanImage[k], i, j, InWidth) == true)//角點偵測
					{
						//將此點加入特徵點資料結構裡面
						AddFeature(kaidanImage[k], i, j, r, insize, InWidth, (k + 1));
					}
				}
			}
		}
	}
}
//*** 檢查是否為極值 ***//
bool SIFT::checkmaxmin(float ** doImage, int nowx, int nowy, int nowz,int InWidth)
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
				if (!(maxc | minc))
					return false;
			}
		}
	}
	value = abs(value);
	if (abs(value) < 0.015)
		return false;
	else
		return true;
}
//*** 高斯模糊 ***//
void SIFT::GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma)
{
	

	float musk[3];
	GusPersent(musk, 3);
	float block[3];
	float total;
	//做橫向的高斯模糊
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (i == 0)
			{
				block[0] = doImage[inz - 1][j*InWidth + i + 0];
				block[1] = doImage[inz - 1][j*InWidth + i + 0];
				block[2] = doImage[inz - 1][j*InWidth + i + 1];
			}
			else if (i == (InWidth - 1))
			{
				block[0] = doImage[inz - 1][j*InWidth + i - 1];
				block[1] = doImage[inz - 1][j*InWidth + i + 0];
				block[2] = doImage[inz - 1][j*InWidth + i + 0];
			}
			else
			{
				block[0] = doImage[inz - 1][j*InWidth + i - 1];
				block[1] = doImage[inz - 1][j*InWidth + i + 0];
				block[2] = doImage[inz - 1][j*InWidth + i + 1];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v]*musk[v];
			}
			doImage[inz][j*InWidth + i] = total;
		}
	}
	//做縱向的高斯模糊
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (j == 0)
			{
				block[0] = doImage[inz][(j + 0)*InWidth + i];
				block[1] = doImage[inz][(j + 0)*InWidth + i];
				block[2] = doImage[inz][(j + 1)*InWidth + i];
			}
			else if (j == (InHeight - 1))
			{
				block[0] = doImage[inz][(j - 1)*InWidth + i];
				block[1] = doImage[inz][(j + 0)*InWidth + i];
				block[2] = doImage[inz][(j + 0)*InWidth + i];
			}
			else
			{
				block[0] = doImage[inz][(j - 1)*InWidth + i];
				block[1] = doImage[inz][(j + 0)*InWidth + i];
				block[2] = doImage[inz][(j + 1)*InWidth + i];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v] * musk[v];
			}
			doImage[inz][j*InWidth + i] = total;
		}
	}
	
}
//*** 計算高斯遮罩的機率(一維) ***//
void SIFT::GusPersent(float *inarray, int line,float sigma)
{
	float *block = new float[line];
	float total = 0;
	int origin = (line + 1) / 2 - 1;
	for (int i = origin; i < line; i++)
	{
		block[i] = pow(E, -(i - origin)*(i - origin) / (2 * sigma*sigma)) / (sqrt(2 * PI*sigma*sigma));
		block[origin*2-i] = pow(E, -(i - origin)*(i - origin) / (2 * sigma*sigma)) / (sqrt(2 * PI*sigma*sigma));
	}
	for (int i = 0; i < line; i++)
	{
		total += block[i];
	}
	for (int i = 0; i < line; i++)
	{
		inarray[i] = block[i] / total;
	}
	delete[] block;
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
//*** 印出箭頭 ***//
void SIFT::display()
{
	this->Height;
	Featureptr pictureS = FeatureStart;
	while (pictureS->nextptr != NULL)
	{
		pictureS = pictureS->nextptr;
		int x, y;
		x = 1;
		while (true)
		{
			y = x * tan(pictureS->sita);
			if ((x*x + y*y) <= (pictureS->mm)*8000.0)
			{
				color[pictureS->y + y][pictureS->x + x].rgbtBlue = 255.0;
				color[pictureS->y + y][pictureS->x + x].rgbtGreen = 0.0;
				color[pictureS->y + y][pictureS->x + x].rgbtRed = 0.0;
			}
			else
			{
				break;
			}
			++x;
		}
	}
}
//*** 在特徵點串列後頭加入新的特徵點 ***//
void SIFT::AddnewFeaturestruct(int Inx, int Iny, float Inm, int Insita)
{
	Featureptr newnode = new Feature;
	newnode->x = Inx;
	newnode->y = Iny;
	newnode->mm = Inm;
	newnode->sita = Insita;
	newnode->nextptr = NULL;
	FeatureEnd->nextptr = newnode;
	FeatureEnd = newnode;
}
//添加特徵點
void SIFT::AddFeature(float * doImage, int Inx, int Iny,int Inr,int Insize,int InWidth,float sigma)
{
	int newLength = Inr * 2 + 1;
	float *musk = new float[newLength];//高斯分布的機率大小
	GusPersent(musk, (Inr * 2 + 1), sigma);
	float direction[36] = {0};//儲存36個方向的加總值
	float *m, *sita;//建構出遮罩大小的m矩陣與sita矩陣
	m = new float[newLength * newLength];
	sita = new float[newLength * newLength];
	int starty, endy;
	int startx, endx;
	starty = Iny - Inr; endy = Iny + Inr;
	startx = Inx - Inr; endx = Inx + Inr;
	for (int j = starty; j <= endy; j++)
	{
		for (int i = startx; i <= endx; i++)
		{
			//做m、sita計算
			float Ldx, Ldy;
			Ldx = doImage[(j + 0)*InWidth + (i + 1)] - doImage[(j + 0)*InWidth + (i - 1)];
			Ldy = doImage[(j + 1)*InWidth + (i + 0)] - doImage[(j - 1)*InWidth + (i + 0)];
			m[(j - starty)*newLength + (i - startx)] = sqrt(Ldx*Ldx + Ldy*Ldy);
			if (Ldx == 0)
			{
				if (Ldy >= 0) sita[(j - starty)*newLength + (i - startx)] = 90.0;
				else sita[(j - starty)*newLength + (i - startx)] = 270.0;
			}
			else
			{
				float getsita = tanh(Ldy / Ldx)*180.0 / PI;
				if (getsita < 0) getsita = 360.0 + getsita;
				sita[(j - starty)*newLength + (i - startx)] = getsita;
			}
		}
	}
	for (int j = 0; j < newLength; j++)
	{
		for (int i = 0; i < newLength; i++)
		{
			//做橫向高斯
			m[j * newLength + i] *= musk[i];
		}
	}
	for (int j = 0; j < newLength; j++)
	{
		for (int i = 0; i < newLength; i++)
		{
			//做縱向高斯
			m[j * newLength + i] *= musk[j];
		}
	}
	//計算矩陣內各方向的加總
	for (int j = 0; j < newLength; j++)
	{
		for (int i = 0; i < newLength; i++)
		{
			int usesita = (sita[j * newLength + i]) / 10;
			direction[usesita] += m[j * newLength + i];
		}
	}
	int sitafoam[36];
	int changesita;
	float changedirection;
	for (int i = 0; i < 36; i++)
	{
		sitafoam[i] = i * 10;
	}
	//將個角度內各大小地值做排序
	for (int j = 0; j < 35; j++)
	{
		for (int i = 0; i < 35-j; i++)
		{
			if(direction[i] <= direction[i+1])
			{
				changesita = sitafoam[i];
				sitafoam[i] = sitafoam[i + 1];
				sitafoam[i + 1] = changesita;

				changedirection = direction[i];
				direction[i] = direction[i + 1];
				direction[i + 1] = changedirection;
			}
		}
	}
	//找出主副方向並加入特徵點結構裡面
	AddnewFeaturestruct(Inx / Insize, Iny / Insize, direction[0], sitafoam[0]);
	int add = 1;
	while (direction[add] >= 0.8*direction[0])
	{
		AddnewFeaturestruct(Inx / Insize, Iny / Insize, direction[add], sitafoam[add]);
		++add;
	}
	delete[] musk;
	delete[] m;
	delete[] sita;
}
//完成一連串尋找特徵點的動作
void SIFT::doFeature()
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
		KaIDaN[0] = new float[inWidth*inHeight];
		ZoomInOut(KaIDaN[0], inWidth, inHeight);
		for (int i = 1; i < 5; ++i)
		{
			KaIDaN[i] = new float[inWidth*inHeight];
			GusB(KaIDaN, i, inWidth, inHeight);
		}
		//高斯差
		for (int i = 0; i < 4; i++)
		{
			GusC(KaIDaN, i, inWidth, inHeight);
		}
		//擷取特徵點
		getFeature(KaIDaN, inWidth, inHeight,size);
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
//*** 轉灰階並且正規化(pixel/255) ***//
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
void SIFT::OutBMP(string outname)
{
	int fix;
	/* 計算每列需略過的 bytes 數 */
	if ((Width * 3) % 4 != 0)
		fix = 4 - ((Width * 3) % 4);
	else
		fix = 0;

	fstream out;
	outname += ".bmp";
	out.open(outname, ios::out | ios::binary);
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
void SIFT::OutRAW(string outname)
{
	fstream out;
	outname += ".raw";
	out.open(outname, ios::out | ios::binary);

	for (int j = 0; j < Height; ++j)
	{
		for (int i = 0; i < Width; ++i)
		{
			out.write((char*)&color[j][i], sizeof(RGBTRIPLE));
		}
	}
}
