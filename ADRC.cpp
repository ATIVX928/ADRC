#include "Adrc.h"
#include "math.h"

//const float ADRC_Unit[17] =

/*TD����΢����   �Ľ�����TD,h0=N*h      ����״̬�۲���ESO           �Ŷ�����     ���������*/
/*  r     h      N                  beta_01   beta_02    beta_03     b0     beta_0  beta_1     beta_2     N1     C    alpha1  alpha2  zeta  b alpha0*/
//{ 10000 ,0.001f , 3,               300,      4000,      10000,     0.001f,    0.00f,   5.0f,      0.010f,    5,    5,    0.75f,   1.25f,    50,    0,-0.2f };


float Constrain_Float(float amt, float low, float high) {
	return ((amt) < (low) ? (low) : ((amt) > (high) ? (high) : (amt)));
}

int Sign_ADRC(float Input)
{
	int output = 0;
	if (Input > 1E-6f) output = 1;
	else if (Input < -1E-6f) output = -1;
	else output = 0;
	return output;
}

int Fsg_ADRC(float x, float d)
{
	int output = 0;
	output = (Sign_ADRC(x + d) - Sign_ADRC(x - d)) / 2.f;
	return output;
}


/*void ADRC_Init(Fhan_Data* fhan_Input1)
{
	fhan_Input1->r = ADRC_Unit[0];
	fhan_Input1->h = ADRC_Unit[1];
	fhan_Input1->N0 = (int)(ADRC_Unit[2]);
	fhan_Input1->beta_01 = ADRC_Unit[3];
	fhan_Input1->beta_02 = ADRC_Unit[4];
	fhan_Input1->beta_03 = ADRC_Unit[5];
	fhan_Input1->b0 = ADRC_Unit[6];
	fhan_Input1->beta_0 = ADRC_Unit[7];
	fhan_Input1->beta_1 = ADRC_Unit[8];
	fhan_Input1->beta_2 = ADRC_Unit[9];
	fhan_Input1->N1 = (int)(ADRC_Unit[10]);
	fhan_Input1->c = ADRC_Unit[11];

	fhan_Input1->alpha1 = ADRC_Unit[12];
	fhan_Input1->alpha2 = ADRC_Unit[13];
	fhan_Input1->zeta = ADRC_Unit[14];
	fhan_Input1->b = ADRC_Unit[15];
	fhan_Input1->alpha0 = ADRC_Unit[16];
	fhan_Input1->e0 = 0;


}*/
void ADRC_Init(Fhan_Data* fhan_Input, float r, float h, float z, float p, float d)
{
	fhan_Input->r = r;
	fhan_Input->h = h;
	fhan_Input->N0 = 3;
	fhan_Input->beta_01 = 1.f/h;
	fhan_Input->beta_02 = 1.f / powf(h, 1.5f) / 1.6f;
	fhan_Input->beta_03 = 1.f / powf(h, 2.2f) / 8.6f;
	fhan_Input->b0 = 0.001f;//�Ŷ���������ʱ������
	fhan_Input->beta_0 = 0;//�����ɻ��ֻ��ڣ�������
	fhan_Input->beta_1 = p;
	fhan_Input->beta_2 = d;

	fhan_Input->alpha1 = 0.8f;
	fhan_Input->alpha2 = 1.5f;
	fhan_Input->zeta = z;
	fhan_Input->b = 0;//ESO�Ŷ���������ʱ������
	fhan_Input->alpha0 = -0.2f;//�����ɻ��ֻ��ڷ����Բ�����������
	fhan_Input->e0 = 0;
}



//ADRC���ٸ���΢����TD���Ľ����㷨fhan
void Fhan_ADRC(Fhan_Data* fhan_Input, float expect_ADRC)//����ADRC���ȹ���
{
	float d = 0, a0 = 0, y = 0, a1 = 0, a2 = 0, a = 0;
	float x1_delta = 0;//ADRC״̬���������
	x1_delta = fhan_Input->x1 - expect_ADRC;//��x1-v(k)���x1�õ���ɢ���¹�ʽ
	fhan_Input->h0 = fhan_Input->N0 * fhan_Input->h;//��h0���h��������ٸ���΢�����ٶȳ�������
	d = fhan_Input->r * fhan_Input->h0 * fhan_Input->h0;//d=rh^2;
	a0 = fhan_Input->h0 * fhan_Input->x2;//a0=h*x2
	y = x1_delta + a0;//y=x1+a0
	a1 = sqrtf(d * (d + 8.f * fabsf(y)));//a1=sqrt(d*(d+8*fabsf(y))])
	a2 = a0 + Sign_ADRC(y) * (a1 - d) / 2;//a2=a0+sign(y)*(a1-d)/2;
	a = (a0 + y) * Fsg_ADRC(y, d) + a2 * (1 - Fsg_ADRC(y, d));
	fhan_Input->fh = -fhan_Input->r * (a / d) * Fsg_ADRC(a, d)
		- fhan_Input->r * Sign_ADRC(a) * (1 - Fsg_ADRC(a, d));//�õ�����΢�ּ��ٶȸ�����
	fhan_Input->x1 += fhan_Input->h * fhan_Input->x2;//�������ٸ���״̬��x1
	fhan_Input->x2 += fhan_Input->h * fhan_Input->fh;//�������ٸ���״̬��΢��x2
}


//ԭ�㸽���������Զε������ݴκ���
float Fal_ADRC(float e, float alpha, float zeta)
{
	//int s = 0;
	float fal_output = 0;
	//s = (Sign_ADRC(e + zeta) - Sign_ADRC(e - zeta)) / 2.f;
	if (fabsf(e) <= zeta)fal_output = e / powf(zeta, 1.f - alpha);
	else fal_output = Sign_ADRC(e) * powf(fabsf(e), alpha);
	//fal_output = e * s / (powf(zeta, 1 - alpha)) + powf(fabsf(e), alpha) * Sign_ADRC(e) * (1 - s);
	return fal_output;
}




/************����״̬�۲���********************/
//״̬�۲�������beta01=1/h  beta02=1/(3*h^2)  beta03=2/(8^2*h^3) ...
void ESO_ADRC(Fhan_Data* fhan_Input)
{
	fhan_Input->e = fhan_Input->z1 - fhan_Input->y;//״̬���

	fhan_Input->fe = Fal_ADRC(fhan_Input->e, 0.5f, fhan_Input->h);//�����Ժ�������ȡ����״̬�뵱ǰ״̬���
	fhan_Input->fe1 = Fal_ADRC(fhan_Input->e, 0.25f, fhan_Input->h);

	/*************��չ״̬������**********/
	fhan_Input->z1 += fhan_Input->h * (fhan_Input->z2 - fhan_Input->beta_01 * fhan_Input->e);
	fhan_Input->z2 += fhan_Input->h * (fhan_Input->z3
		- fhan_Input->beta_02 * fhan_Input->fe
		+ fhan_Input->b * fhan_Input->u);
	//ESO����״̬���ٶ��źţ������Ŷ���������ͳMEMS������Ư�ƽϴ󣬹��ƻ����Ư��
	fhan_Input->z3 += fhan_Input->h * (-fhan_Input->beta_03 * fhan_Input->fe1);
}

void Nolinear_Conbination_ADRC(Fhan_Data* fhan_Input)
{
	float temp_e2 = 0;
	temp_e2 = Constrain_Float(fhan_Input->e2, -300, 300);

	fhan_Input->u0 = fhan_Input->beta_1 * Fal_ADRC(fhan_Input->e1, fhan_Input->alpha1, fhan_Input->zeta)
		//+ fhan_Input->beta_0 * Fal_ADRC(fhan_Input->e0, fhan_Input->alpha0, fhan_Input->zeta)
		+ fhan_Input->beta_2 * Fal_ADRC(temp_e2, fhan_Input->alpha2, fhan_Input->zeta);
}


void ADRC_Control(Fhan_Data* fhan_Input, float expect_ADRC, float feedback_ADRC)
{
	/*�Կ��ſ�������1��*/
	/********
		**
		**
		**
		**
		**
	 ********/
	 /*****
	 ���Ź��ȹ��̣�����Ϊ����������
	 ��TD����΢�����õ���
	 ���������ź�x1����������΢���ź�x2
	 ******/
	Fhan_ADRC(fhan_Input, expect_ADRC);

	/*�Կ��ſ�������2��*/
	/********
			*
			*
	   ****
	 *
	 *
	 ********/
	 /************ϵͳ���ֵΪ��������״̬������ESO����״̬�۲���������*********/
	fhan_Input->y = feedback_ADRC;
	/*****
	����״̬�۲������õ������źŵ�����״̬��
	1��״̬�ź�z1��
	2��״̬�ٶ��ź�z2��
	3��״̬���ٶ��ź�z3��
	����z1��z2������Ϊ״̬������TD΢�ָ������õ���x1,x2�����
	���������Ժ���ӳ�䣬����betaϵ����
	��ϵõ�δ����״̬���ٶȹ����Ŷ�������ԭʼ������u
	*********/
	ESO_ADRC(fhan_Input);//�ͳɱ�MEMS�����Ư�ƣ���չ������z3�����Ư�ƣ�Ŀǰ��ʱδ�뵽�취�����δ�õ�z3
  /*�Կ��ſ�������3��*/
  /********
		 **
	   **
	 **
	   **
		 **
   ********/
   /********״̬������***/
	fhan_Input->e0 += fhan_Input->e1 * fhan_Input->h;//״̬������
	fhan_Input->e1 = fhan_Input->x1 - fhan_Input->z1;//״̬ƫ����
	fhan_Input->e2 = fhan_Input->x2 - fhan_Input->z2;//״̬΢���
	/********�������*******/
   /*
	fhan_Input->u0=//fhan_Input->beta_0*fhan_Input->e0
				  +fhan_Input->beta_1*fhan_Input->e1
				  +fhan_Input->beta_2*fhan_Input->e2;
   */
	Nolinear_Conbination_ADRC(fhan_Input);
	/**********�Ŷ�����*******/
	//fhan_Input->u=fhan_Input->u0
	//             -fhan_Input->z3/fhan_Input->b0;
	//����MEMS������Ư�ƱȽ����أ���beta_03ȡֵ�Ƚϴ�ʱ����ʱ��z3Ư�ƱȽϴ�Ŀǰ�������Ŷ�����������
	fhan_Input->u = Constrain_Float(fhan_Input->u0, -500, 500);
}

