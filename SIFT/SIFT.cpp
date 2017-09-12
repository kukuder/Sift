#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <Windows.h>
#include <string>
#include "SIFT.h"
#define PI 3.14159265358979323846
#define E  2.718281828459045235
#define SIGMA 1.6
/* test */
using namespace std;

struct RGBDATA {
	BYTE rgbtRed;
	BYTE rgbtGreen;
	BYTE rgbtBlue;
};

struct Feature
{
	int x, y;//�U�Ҧb���h���y��
	float mm;//�j��
	int sita;//�]�t�D��V�P�t��V������
	float size;//��
	int kai;//�h
	float sigmaOCT;
	vector <vector <vector <float>>> descrip;
	Feature *nextptr = nullptr;
};
typedef struct Feature* Featureptr;

SIFT::SIFT(char* inname)
{
	Readbmp(inname);
	RGB2Gray();
	FeatureStart = new Feature;
	FeatureEnd = FeatureStart;
	doFeature();
}

SIFT::~SIFT()
{
	for (int i = 0; i < Height; ++i)
	{
		delete[] color[i];
	}
	delete[] color;
	Destroy();
}

void SIFT::Readbmp(char * inname)
{
	int fix;
	/* �}���ɮ� */
	fstream fp;
	fp.open(inname, ios::in | ios::binary);
	/* ��Ū�����Y��T */
	fp.read((char*)&FileHeader, sizeof(BITMAPFILEHEADER));
	fp.read((char*)&InfoHeader, sizeof(BITMAPINFOHEADER));
	cout << FileHeader.bfSize << endl;
	cout << InfoHeader.biSizeImage << endl;
	/* ���o�ϼe�ιϰ� */
	Width = InfoHeader.biWidth;
	Height = InfoHeader.biHeight;
	
	color = new RGBTRIPLE*[Height];
	for (int i = 0; i < Height; ++i)
	{
		color[i] = new RGBTRIPLE[Width];
	}

	cout << "Width : " << Width << endl;
	cout << "Height : " << Height << endl;

	/* �p��C�C�ݲ��L�� bytes �� */
	if ((Width * 3) % 4 != 0)
		fix = 4 - ((Width * 3) % 4);
	else
		fix = 0;

	/* �ϤJ�U������T */
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
		/* ���L�U�C�h�l����T */
		BYTE   ByteBuf;
		for (int n = 0; n<fix; ++n) {
			fp.read((char*)&ByteBuf, sizeof(BYTE));
		}
	}

	//****************************��X
	//OutBMP("out");
	//**************************
}
//*** ���I���� ***//
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
	if ((Tr*Tr / Det) < 12.1)//12.1 �O (r+1)/r  r��10 �ұo����
	{
		return true;
	}
	else
	{
		return false;
	}
}
//*** �^���S�x�I ***//
void SIFT::getFeature(float** kaidanImage, float** DogImage, int InWidth, int InHeight, int insize, float sigma)
{
	for (int k = 1; k <= 2; k++)
	{
		int r = insize * 3 * 1.5*(k + 1);
		for (int j = r + 1; j < (InHeight - r - 1); j++)
		{
			for (int i = r + 1; i < (InWidth - r - 1); i++)
			{
				bool rof = checkmaxmin(DogImage, i, j, k,InWidth);//�����˴�
				if (rof == true)
				{
					if (Hessian(DogImage[k], i, j, InWidth) == true)//���I����
					{
						//�N���I�[�J�S�x�I��Ƶ��c�̭�
						AddFeature(kaidanImage[k], i, j, r, insize, InWidth, pow(sqrt(2), k)*sigma,k);
					}
				}
			}
		}
	}
}
//*** �ˬd�O�_������ ***//
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
//*** �����ҽk ***//
void SIFT::GusB(float** doImage,int inz, int InWidth, int InHeight, float sigma)
{
	float musk[3];
	GusPersent(musk, 3, sigma);
	float block[3];
	float total;
	//����V�������ҽk
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
	//���a�V�������ҽk
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
//*** �p�Ⱚ���B�n�����v(�@��) ***//
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
//*** �����t ***//
void SIFT::GusC(float** FinalImage, float** DogImage, int inz, int InWidth, int InHeight)
{
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			float value1 = FinalImage[inz][j*InWidth + i];
			float value2 = FinalImage[inz+1][j*InWidth + i];
			value1 -= value2;
			DogImage[inz][j*InWidth + i] = value1;
		}
	}
}
//*** �L�X�b�Y ***//
void SIFT::display()
{
	int mag = 10000.0;//�b�Y���v
	Featureptr pictureS = FeatureStart;
	while (pictureS->nextptr != NULL)
	{
		pictureS = pictureS->nextptr;
		int x, y;
		float roominout = pictureS->size;
		x = 1;
		//���X�b�Y�����u
		while (true)
		{
			if (pictureS->sita == 90 || pictureS->sita == 270) y = 0;
			else y = x * tan((pictureS->sita)*PI / 180.0);

			if ((x*x + y*y) <= (pictureS->mm)*mag)
			{
				int dx, dy;
				if (pictureS->sita > 90 && pictureS->sita <= 270)
				{
					dx = (int)(pictureS->x / roominout) - x;
					dy = (int)(pictureS->y / roominout) - y;
				}
				else
				{
					dx = (int)(pictureS->x / roominout) + x;
					dy = (int)(pictureS->y / roominout) + y;
				}

				if (dx >= 0 && dx < Width && dy >= 0 && dy < Height)
				{
					color[dy][dx].rgbtBlue = 255.0;
					color[dy][dx].rgbtGreen = 0.0;
					color[dy][dx].rgbtRed = 0.0;
				}
			}
			else
			{
				break;
			}
			++x;
		}
		if (pictureS->sita > 90 && pictureS->sita <= 270)
		{
			x = -x;
			y = -y;
		}
		//�e�X�b�Y������׽u
		int xR, yR;
		int xL, yL;
		xR = xL = 1;
		int sitaR = pictureS->sita - 150;
		int sitaL = pictureS->sita + 150;
		if (sitaR < 0) sitaR = 360 + sitaR;
		if (sitaL >= 360) sitaL = sitaL - 360;
		while (true)//�k�U�u
		{
			if (sitaR == 90 || sitaR == 270) yR = 0;
			else yR = xR * tan((sitaR)*PI / 180.0);

			if ((xR*xR + yR*yR) <= (pictureS->mm)*mag / 4.0)
			{
				int dx, dy;
				if (sitaR > 90 && sitaR <= 270)
				{
					dx = pictureS->x / roominout + x - xR;
					dy = pictureS->y / roominout + y - yR;
				}
				else
				{
					dx = pictureS->x / roominout + x + xR;
					dy = pictureS->y / roominout + y + yR;
				}

				if (dx >= 0 && dx < Width && dy >= 0 && dy < Height)
				{
					color[dy][dx].rgbtBlue = 255.0;
					color[dy][dx].rgbtGreen = 0.0;
					color[dy][dx].rgbtRed = 0.0;
				}
			}
			else
			{
				break;
			}
			++xR;
		}

		while (true)//���U�u
		{
			if (sitaL == 90 || sitaL == 270) yL = 0;
			else yL = xL * tan((sitaL)*PI / 180.0);

			if ((xL*xL + yL*yL) <= (pictureS->mm)*mag / 4.0)
			{
				int dx, dy;
				if (sitaL > 90 && sitaL <= 270)
				{
					dx = pictureS->x / roominout + x - xL;
					dy = pictureS->y / roominout + y - yL;
				}
				else
				{
					dx = pictureS->x / roominout + x + xL;
					dy = pictureS->y / roominout + y + yL;
				}

				if (dx >= 0 && dx < Width && dy >= 0 && dy < Height)
				{
					color[dy][dx].rgbtBlue = 255.0;
					color[dy][dx].rgbtGreen = 0.0;
					color[dy][dx].rgbtRed = 0.0;
				}
			}
			else
			{
				break;
			}
			++xL;
		}
	}
}
//*** �b�S�x�I��C���Y�[�J�s���S�x�I ***//
void SIFT::AddnewFeaturestruct(int Inx, int Iny, int Insize, int kai, int sigmaOCT, float Inm, int Insita)
{
	Featureptr newnode = new Feature;
	newnode->x = Inx;
	newnode->y = Iny;
	newnode->mm = Inm;
	newnode->sita = Insita;
	newnode->size = Insize;
	newnode->kai = kai;
	newnode->sigmaOCT = sigmaOCT;
	newnode->nextptr = NULL;
	FeatureEnd->nextptr = newnode;
	FeatureEnd = newnode;
}
//*** �K�[�S�x�I ***//
void SIFT::AddFeature(float * doImage, int Inx, int Iny,int Inr,int Insize,int InWidth,float sigma, int kai)
{
	int newLength = Inr * 2 + 1;
	float *musk = new float[newLength];//�������������v�j�p
	GusPersent(musk, (Inr * 2 + 1), sigma);
	float direction[36] = {0};//�x�s36�Ӥ�V���[�`��
	float *mm, *sita;//�غc�X�B�n�j�p��m�x�}�Psita�x�}
	mm = new float[newLength * newLength];
	sita = new float[newLength * newLength];
	int starty, endy;
	int startx, endx;
	starty = Iny - Inr; endy = Iny + Inr;
	startx = Inx - Inr; endx = Inx + Inr;
	for (int j = starty; j <= endy; j++)
	{
		for (int i = startx; i <= endx; i++)
		{
			//��m�Bsita�p��
			mm[(j - starty)*newLength + (i - startx)] = getM(doImage[(j + 0)*InWidth + (i - 1)], doImage[(j + 0)*InWidth + (i + 1)], doImage[(j - 1)*InWidth + (i + 0)], doImage[(j + 1)*InWidth + (i + 0)]);
			sita[(j - starty)*newLength + (i - startx)] = getSita(doImage[(j + 0)*InWidth + (i - 1)], doImage[(j + 0)*InWidth + (i + 1)], doImage[(j - 1)*InWidth + (i + 0)], doImage[(j + 1)*InWidth + (i + 0)]);
			sita[(j - starty)*newLength + (i - startx)] = sita[(j - starty)*newLength + (i - startx)] * 180.0 / PI;
		}
	}
	for (int j = 0; j < newLength; j++)
	{
		for (int i = 0; i < newLength; i++)
		{
			//����V����
			mm[j * newLength + i] *= musk[i];
		}
	}
	for (int j = 0; j < newLength; j++)
	{
		for (int i = 0; i < newLength; i++)
		{
			//���a�V����
			mm[j * newLength + i] *= musk[j];
		}
	}
	//�p��x�}���U��V���[�`
	for (int j = 0; j < newLength; j++)
	{
		for (int i = 0; i < newLength; i++)
		{
			int usesita = (sita[j * newLength + i]) / 10;
			direction[usesita] += mm[j * newLength + i];
		}
	}
	int sitafoam[36];
	int changesita;
	float changedirection;
	for (int i = 0; i < 36; i++)
	{
		sitafoam[i] = i * 10;
	}
	//�N�Ө��פ��U�j�p�a�Ȱ��Ƨ�
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
	//��X�D�Ƥ�V�å[�J�S�x�I���c�̭�
	AddnewFeaturestruct(Inx, Iny, Insize, kai, sigma, direction[0], sitafoam[0]);
	int add = 1;
	while (direction[add] >= 0.8*direction[0])
	{
		AddnewFeaturestruct(Inx, Iny, Insize, kai, sigma, direction[add], sitafoam[add]);
		++add;
	}
	delete[] musk;
	delete[] mm;
	delete[] sita;
}
//*** �����@�s��M��S�x�I���ʧ@ ***//
void SIFT::doFeature()
{
	//************** �H�Ϥ��j�p�p��X�����r�𰵦h�� ****************//
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
	//********************* ��j�Y�p�H�ΰ����ҽk *************************//
	float* FirstPicture;
	int nowkai = 0;
	FirstPicture = new float[Width*Height];
	GetGrayPicture(FirstPicture);
	for (float size = 1; size >= 1; size/=2)
	{
		float** KaIDaN;
		float** KaIDaNGray;
		KaIDaN = new float*[5];
		int inWidth = Width*size;
		int inHeight = Height*size;
		//�����ҽk
		KaIDaN[0] = new float[inWidth*inHeight];
		ZoomInOut(KaIDaN[0], inWidth, inHeight);
		//���տ�X
		OutBMPCheck("�Ƕ���", inWidth, inHeight, KaIDaN);//���տ�X
		if (nowkai == 0)
		{
			Gusto(KaIDaN[0], inWidth, inHeight, SIGMA);
		}
		else
		{
			Gusto(KaIDaN[0], inWidth, inHeight, SIGMA*nowkai * 2);
		}
		for (int i = 1; i < 5; ++i)
		{
			KaIDaN[i] = new float[inWidth*inHeight];
			GusB(KaIDaN, i, inWidth, inHeight, SIGMA*pow(sqrt(2),i));
		}
		//���տ�X
		OutRAW("out");
		OutBMPCheck("����", inWidth, inHeight, KaIDaN);//���տ�X
		//�����t
		KaIDaNGray = new float*[3];
		for (int i = 0; i < 4; i++)
		{
			KaIDaNGray[i] = new float[inWidth*inHeight];
		}
		for (int i = 0; i < 4; i++)
		{
			GusC(KaIDaN, KaIDaNGray, i, inWidth, inHeight);
		}
		//���տ�X
		OutBMPCheck("�t����", inWidth, inHeight, KaIDaNGray);//���տ�X
		//�^���S�x�I
		FeatureNow = FeatureEnd;
		getFeature(KaIDaN, KaIDaNGray, inWidth, inHeight, size, SIGMA);
		//���ͯS�x�y�z�l
		FeatureDescrip(KaIDaN);
		//����O����
		for (int i = 0; i < 5; i++)
		{
			delete[] KaIDaN[i];
		}
		for (int i = 0; i < 4; i++)
		{
			delete[] KaIDaNGray[i];
		}
	}
	delete[] FirstPicture;
}
//*** ���ͯS�x�y�z�l ***//
void SIFT::FeatureDescrip(float** kaidaImag)
{
	//�s�W�@��4*4*8���S�x�y�z�Ŷ�
	vector <vector <vector <float>>> descripgroup;//�x�s�S�x�y�z�l
	descripgroup.resize(4);
	for (int v = 0; v < 4; v++)
	{
		descripgroup[v].resize(4);
		for (int j = 0; j < 4; j++)
		{
			descripgroup[v][j].resize(8);
		}
	}
	for (int v = 0; v < 4; v++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < 8; i++)
			{
				descripgroup[v][j][i] = 0;
			}
		}
	}
	float newXYArray[16 * 16] = { 0 };
	Featureptr FeatureS = FeatureNow;
	while (FeatureS->nextptr != nullptr)
	{
		FeatureS = FeatureS->nextptr;
		//���
		int radius = (FeatureS->sigmaOCT * 3 * sqrt(2)*(5) + 1) / 2;
		int inHeight = Height*FeatureS->size;
		int inWidth = Width*FeatureS->size;
		for (int j = FeatureS->y - radius ; j < FeatureS->y + radius; j++)
		{
			for (int i = FeatureS->x - radius; i < FeatureS->x + radius; i++)
			{
				if (j < 1 || j >= inHeight - 1 || i < 1 || i >= inWidth - 1)
				{
					continue;
				}
				int newx, newy;
				newx = i - FeatureS->x;
				newy = j - FeatureS->y;
				CoordinateChange(&newx, &newy, (FeatureS->sita * PI / 180.0));
				float blockx = newx / ((radius*sqrt(2) / 2.0)) * 2.0 + 3.0;// ���d�򸨦b1~4������
				float blocky = newy / ((radius*sqrt(2) / 2.0)) * 2.0 + 3.0;// ���d�򸨦b1~4������
				float sita, mm;
				mm = getM(kaidaImag[FeatureS->kai][j*inWidth + i - 1], kaidaImag[FeatureS->kai][j*inWidth + i + 1], kaidaImag[FeatureS->kai][(j - 1)*inWidth + i], kaidaImag[FeatureS->kai][(j + 1)*inWidth + i]);
				sita = getSita(kaidaImag[FeatureS->kai][j*inWidth + i - 1], kaidaImag[FeatureS->kai][j*inWidth + i + 1], kaidaImag[FeatureS->kai][(j - 1)*inWidth + i], kaidaImag[FeatureS->kai][(j + 1)*inWidth + i]);
				float unit = 2 * PI / 8;

				float weix0, weix1;
				float weiy0, weiy1;
				weix1 = blockx - (int)blockx; weix0 = 1 - weix1;
				weiy1 = blocky - (int)blocky; weiy0 = 1 - weiy1;

				//****** �Hfirstorder�N�Ȥ��O�s�Jdescripgroup�e���̭�
				int bx0, bx1;
				int by0, by1;
				bx0 = blockx;
				bx1 = bx0 + 1;
				by0 = blocky;
				by1 = by0 + 1;

				float RU, RD, LU, LD;
				float SitaGroup;
				
				SitaGroup = sita / unit;

				if (bx0 >= 1 && bx0 <= 4 && by0 >= 1 && by0 <= 4)//RU
				{
					RU = weix0 * weiy0 * mm;
					descripgroup[by0 - 1][bx0 - 1][SitaGroup] += RU;
				}
				if (bx0 >= 1 && bx0 <= 4 && by1 >= 1 && by1 <= 4)//RD
				{
					RD = weix0 * weiy1 * mm;
					descripgroup[by1 - 1][bx0 - 1][SitaGroup] += RD;
				}
				if (bx1 >= 1 && bx1 <= 4 && by0 >= 1 && by0 <= 4)//LU
				{
					LU = weix1 * weiy0 * mm;
					descripgroup[by0 - 1][bx1 - 1][SitaGroup] += LU;
				}
				if (bx1 >= 1 && bx1 <= 4 && by1 >= 1 && by1 <= 4)//LD
				{
					LD = weix1 * weiy1 * mm;
					descripgroup[by1 - 1][bx1 - 1][SitaGroup] += LD;
				}

			}
		}
		FeatureS->descrip = descripgroup;
	}
}
//*** ��j�Y�p ***//
void SIFT::ZoomInOut(float* doImage, int InWidth, int InHeight)
{
	float XX1, XX2, XXYY;
	float dx, dy;
	float Io, Jo;
	for (int j = 0; j < InHeight - 1; ++j)
	{
		Jo = j * (Height - 1) / (InHeight - 1);
		dy = 1 - ((j * (Height - 1) / (InHeight - 1)) - Jo);
		for (int i = 0; i < InWidth - 1; ++i)
		{
			Io = i * (Width - 1) / (InWidth - 1);
			dx = 1 - ((i * (Width - 1) / (InWidth - 1)) - Io);
			XX1 = gray[(int)((Io + 0) + (Jo + 0)*Width)] * dx + gray[(int)((Io + 1) + (Jo + 0)*Width)] * (1 - dx);
			XX2 = gray[(int)((Io + 0) + (Jo + 1)*Width)] * dx + gray[(int)((Io + 1) + (Jo + 1)*Width)] * (1 - dx);
			XXYY = XX1*dy + XX2*(1 - dy);
			doImage[j*InWidth + i] = XXYY;
		}
		XX1 = gray[(int)((Jo + 1)*Width - 1)];
		XX2 = gray[(int)((Jo + 2)*Width - 1)];
		XXYY = XX1*dy + XX2*(1 - dy);
		doImage[(j + 1)*InWidth - 1] = XXYY;
	}
	for (int i = 0; i < InWidth - 1; i++)
	{
		Io = i * (Width - 1) / (InWidth - 1);
		dx = 1 - ((i * (Width - 1) / (InWidth - 1)) - Io);
		XX1 = gray[(int)((Height - 1)*Width + (Io + 0))];
		XX2 = gray[(int)((Height - 1)*Width + (Io + 1))];
		XXYY = XX1*dx + XX2*(1 - dx);
		doImage[(InHeight - 1)*InWidth + i] = XXYY;
	}
	doImage[InHeight*InWidth - 1] = gray[Height*Width - 1];
}
//*** ��Ƕ��åB���W��(pixel/255) ***//
void SIFT::RGB2Gray()
{
	//gray.resize(Height*Width);
	gray = new float[Height * Width];
	for (int j = 0; j < Height; ++j)
	{
		for (int i = 0; i < Width; ++i)
		{
			gray[j*Width + i] = color[j][i].rgbtRed*0.299 + color[j][i].rgbtGreen*0.587 + color[j][i].rgbtBlue*0.114;
			//gray[j*Width + i] /= 255.0;
		}
	}
}
//*** ���o���W�ƪ��Ƕ��Ϥ� ***//
void SIFT::GetGrayPicture(float* doImage)
{
	for (int i = 0; i < Height*Width; i++)
	{
		doImage[i] = 0;
		doImage[i] = gray[i];
	}
}
//*** �N�ɮ׿�X��BMP�榡 ***//
void SIFT::OutBMP(string outname)
{
	int fix;
	/* �p��C�C�ݲ��L�� bytes �� */
	int check = (Width * 3) % 4;
	if (check != 0)
		fix = 4 - check;
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
		/* ���L�U�C�h�l����T */
		for (int n = 0; n<fix; ++n) {
			BYTE   ByteBuf;
			out.write((char*)&ByteBuf, sizeof(ByteBuf));
		}
	}
}
//*** �N�ɮ׿�X��RAW�榡 ***//
void SIFT::OutRAW(string outname)
{
	fstream out;
	outname += ".raw";
	out.open(outname, ios::out | ios::binary);

	for (int j = 0; j < Height * Width; ++j)
	{
		float value = gray[j];
		//value *= 255.0;
		if (value > 255)value = 255.0;
		else if (value < 0) value = 0.0;
		out.write((char*)&value, sizeof(float));
	}
}
//*** �쵲��C�Ѻc ***//
void SIFT::Destroy()
{
	Featureptr nowptr = FeatureStart;
	while (nowptr == nullptr)
	{
		nowptr = FeatureStart->nextptr;
		delete[] FeatureStart;
		FeatureStart = nowptr;
	}
}
//*** �����ҽk(�q�Ϊ�) ***//
void SIFT::Gusto(float * doImage, int InWidth, int InHeight, float sigma)
{
	float musk[3];
	GusPersent(musk, 3, sigma);
	float block[3];
	float total;
	//����V�������ҽk
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (i == 0)
			{
				block[0] = doImage[j*InWidth + i + 0];
				block[1] = doImage[j*InWidth + i + 0];
				block[2] = doImage[j*InWidth + i + 1];
			}
			else if (i == (InWidth - 1))
			{
				block[0] = doImage[j*InWidth + i - 1];
				block[1] = doImage[j*InWidth + i + 0];
				block[2] = doImage[j*InWidth + i + 0];
			}
			else
			{
				block[0] = doImage[j*InWidth + i - 1];
				block[1] = doImage[j*InWidth + i + 0];
				block[2] = doImage[j*InWidth + i + 1];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v] * musk[v];
			}
			doImage[j*InWidth + i] = total;
		}
	}
	//���a�V�������ҽk
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (j == 0)
			{
				block[0] = doImage[(j + 0)*InWidth + i];
				block[1] = doImage[(j + 0)*InWidth + i];
				block[2] = doImage[(j + 1)*InWidth + i];
			}
			else if (j == (InHeight - 1))
			{
				block[0] = doImage[(j - 1)*InWidth + i];
				block[1] = doImage[(j + 0)*InWidth + i];
				block[2] = doImage[(j + 0)*InWidth + i];
			}
			else
			{
				block[0] = doImage[(j - 1)*InWidth + i];
				block[1] = doImage[(j + 0)*InWidth + i];
				block[2] = doImage[(j + 1)*InWidth + i];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v] * musk[v];
			}
			doImage[j*InWidth + i] = total;
		}
	}

}
//*** ��XBMP(�ˬd��) ***//
void SIFT::OutBMPCheck(std::string outname , int inwidth, int inheight, float ** doImage)
{
	BITMAPFILEHEADER checkFileHeader = FileHeader;
	BITMAPINFOHEADER checkInfoHeader = InfoHeader;

	checkInfoHeader.biSizeImage = inwidth * inheight;
	checkInfoHeader.biWidth = inwidth;
	checkInfoHeader.biHeight = inheight;
	checkFileHeader.bfSize = checkInfoHeader.biSizeImage + 54;

	int fix;
	/* �p��C�C�ݲ��L�� bytes �� */
	int check = (inwidth * 3) % 4;
	if (check != 0)
		fix = 4 - check;
	else
		fix = 0;

	fstream out;
	outname += ".bmp";
	out.open(outname, ios::out | ios::binary);
	out.write((char*)&checkFileHeader, sizeof(BITMAPFILEHEADER));
	out.write((char*)&checkInfoHeader, sizeof(BITMAPINFOHEADER));

	for (int j = inheight - 1; j >= 0; --j)
	{
		for (int i = 0; i < inwidth; ++i)
		{
			RGBTRIPLE rgb;
			rgb.rgbtBlue = doImage[0][i + j * inwidth];
			if (rgb.rgbtBlue < 0.0) rgb.rgbtBlue = 0.0;
			else rgb.rgbtBlue *= 255.0;

			if (rgb.rgbtBlue > 255.0)rgb.rgbtBlue = 255.0;

			rgb.rgbtGreen = rgb.rgbtRed = rgb.rgbtBlue;

			out.write((char*)&rgb, sizeof(RGBTRIPLE));
		}
		/* ���L�U�C�h�l����T */
		for (int n = 0; n<fix; ++n) {
			BYTE   ByteBuf;
			out.write((char*)&ByteBuf, sizeof(ByteBuf));
		}
	}
}


//*** �y���ഫ ***//
void CoordinateChange(int* deltaX, int* deltaY, float sita)
{
	int newxx, newyy;
	newxx = cos(sita)*(*deltaX) - sin(sita)*(*deltaY);
	newyy = sin(sita)*(*deltaX) + cos(sita)*(*deltaY);
	(*deltaX) = newxx;
	(*deltaY) = newyy;
}
//*** �Dm(�j��)�� ***//
float getM(float x0, float x1, float y0, float y1)
{
	float mm;
	mm = (x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0);
	mm = sqrt(mm);
	return mm;
}
//*** �Dsita��(�|��) ***//
float getSita(float x0, float x1, float y0, float y1)
{
	float sita;
	if (x0 == x1)
	{
		if (y1 > y0) sita = PI / 2;
		else sita = PI * 3 / 2;
	}
	else
	{
		sita = tanh((y1 - y0) / (x1 - x0));
		if (x1 < x0) sita += PI;
		if (sita < 0) sita = 2 * PI + sita;
		if (sita >= 2 * PI) sita -= 2 * PI;
	}

	

	return sita;
}







Stitching::Stitching(BITMAPFILEHEADER Filedata, BITMAPINFOHEADER Infodata)
{
	FileHeader = Filedata;
	InfoHeader = Infodata;

	
	InfoHeader.biWidth = Infodata.biWidth * 2;
	Width = Infodata.biWidth * 2;
	Height = InfoHeader.biHeight;

	int widthsize;
	if (Width % 4 != 0)
	{
		widthsize = Width * 3 + (4 - Width % 4);
	}

	InfoHeader.biSizeImage = widthsize * Height;
	FileHeader.bfSize = InfoHeader.biSizeImage + 54;
}

void Stitching::OutBMP(string outname)
{
	int fix;
	/* �p��C�C�ݲ��L�� bytes �� */
	int check = (Width * 3) % 4;
	if (check != 0)
		fix = 4 - check;
	else
		fix = 0;

	fstream out;
	outname += ".bmp";
	out.open(outname, ios::out | ios::binary);
	out.write((char*)&FileHeader, sizeof(BITMAPFILEHEADER));
	out.write((char*)&InfoHeader, sizeof(BITMAPINFOHEADER));

	for (int j = Height - 1; j >= 0; --j)
	{
		for (int i = 0; i < Width; ++i)
		{
			RGBDATA rgb;
			rgb.rgbtBlue = color[j][i].rgbtBlue;
			rgb.rgbtGreen = color[j][i].rgbtGreen;
			rgb.rgbtRed = color[j][i].rgbtRed;
			out.write((char*)&rgb, sizeof(RGBTRIPLE));
		}
		/* ���L�U�C�h�l����T */
		for (int n = 0; n<fix; ++n) {
			BYTE   ByteBuf;
			out.write((char*)&ByteBuf, sizeof(ByteBuf));
		}
	}
}
//*** �N�ɮ׿�X��RAW�榡 ***//
void Stitching::OutRAW(string outname)
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
//*** 
void Stitching::WriteImage(RGBTRIPLE** image1, RGBTRIPLE** image2)
{
	color = new RGBTRIPLE*[Height];
	for (int i = 0; i < Height; ++i)
	{
		color[i] = new RGBTRIPLE[Width];
	}

	for (int j = 0; j < Height; j++)
	{
		// �g�J�Ĥ@�i�Ϫ���T
		for (int i = 0; i < Width/2; i++)
		{
			color[j][i] = image1[j][i];
		}
		// �g�J�ĤG�i�Ϫ���T
		for (int i = Width / 2; i < Width; i++)
		{
			color[j][i] = image2[j][i - (Width / 2)];
		}
	}
}
