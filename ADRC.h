#ifndef _ADRC_H_
#define _ADRC_H_

typedef struct ADRC_Parameter
{
	/*****���Ź��ȹ���*******/
	float x1=0;//����΢����״̬��
	float x2=0;//����΢����״̬��΢����
	float r=0;//ʱ��߶�
	float h=0;//ADRCϵͳ����ʱ��
	int N0=0;//����΢��������ٶȳ���h0=N*h

	float h0=0;
	float fh=0;//����΢�ּ��ٶȸ�����
	/*****����״̬�۲���*******/
	/******��ϵͳ���y������u�����ٹ���ϵͳ״̬���Ŷ�*****/
	float z1=0;
	float z2=0;
	float z3=0;//���ݿ��ƶ����������������ȡ���Ŷ���Ϣ
	float e=0;//ϵͳ״̬���
	float y=0;//ϵͳ�����
	float fe=0;
	float fe1=0;
	float beta_01=0;
	float beta_02=0;
	float beta_03=0;
	float b=0;


	/**********ϵͳ״̬������*********/
	float e0=0;//״̬��������
	float e1=0;//״̬ƫ��
	float e2=0;//״̬��΢����
	float u0=0;//���������ϵͳ���
	float u=0;//���Ŷ�����������
	float b0=0;//�Ŷ�����

	/*********��һ�������ʽ*********/
	float beta_0=0;//����
	float beta_1=0;//��������ϲ���
	float beta_2=0;//u0=beta_1*e1+beta_2*e2+(beta_0*e0);
	/*********�ڶ��������ʽ*********/
	float alpha0=0;
	float alpha1=0;//u0=beta_1*fal(e1,alpha1,zeta)+beta_2*fal(e2,alpha2,zeta)
	float alpha2=0;//0<alpha1<1<alpha2
	float zeta=0;//���Զε����䳤��
	/*********�����������ʽ*********/
	//float h1;//u0=-fhan(e1,e2,r,h1);
	//int N1;//����΢��������ٶȳ���h0=N*h
	/*********�����������ʽ*********/
	//float c;//u0=-fhan(e1,c*e2*e2,r,h1);


}Fhan_Data;


//��������ADRC_Init
//���ã���ʼ��
//����r��΢�ָ�����������ԽС�����ٶ�Խ��������Խƽ��
//����h��ʱ�䲽��
//����z������������Ա߽�
//����p��������PID��������
//����d��������PID΢�ֻ���
void ADRC_Init(Fhan_Data* fhan_Input, float r, float h, float z, float p, float d);
void Fhan_ADRC(Fhan_Data* fhan_Input, float expect_ADRC);
void ADRC_Control(Fhan_Data* fhan_Input, float expect_ADRC, float feedback);

extern Fhan_Data ADRC_dx, ADRC_dy, ADRC_dz;
extern Fhan_Data dz_td;
#endif

