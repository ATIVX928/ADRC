#ifndef _ADRC_H_
#define _ADRC_H_

typedef struct ADRC_Parameter
{
	/*****安排过度过程*******/
	float x1=0;//跟踪微分期状态量
	float x2=0;//跟踪微分期状态量微分项
	float r=0;//时间尺度
	float h=0;//ADRC系统积分时间
	int N0=0;//跟踪微分器解决速度超调h0=N*h

	float h0=0;
	float fh=0;//最速微分加速度跟踪量
	/*****扩张状态观测器*******/
	/******已系统输出y和输入u来跟踪估计系统状态和扰动*****/
	float z1=0;
	float z2=0;
	float z3=0;//根据控制对象输入与输出，提取的扰动信息
	float e=0;//系统状态误差
	float y=0;//系统输出量
	float fe=0;
	float fe1=0;
	float beta_01=0;
	float beta_02=0;
	float beta_03=0;
	float b=0;


	/**********系统状态误差反馈率*********/
	float e0=0;//状态误差积分项
	float e1=0;//状态偏差
	float e2=0;//状态量微分项
	float u0=0;//非线性组合系统输出
	float u=0;//带扰动补偿后的输出
	float b0=0;//扰动补偿

	/*********第一种组合形式*********/
	float beta_0=0;//线性
	float beta_1=0;//非线性组合参数
	float beta_2=0;//u0=beta_1*e1+beta_2*e2+(beta_0*e0);
	/*********第二种组合形式*********/
	float alpha0=0;
	float alpha1=0;//u0=beta_1*fal(e1,alpha1,zeta)+beta_2*fal(e2,alpha2,zeta)
	float alpha2=0;//0<alpha1<1<alpha2
	float zeta=0;//线性段的区间长度
	/*********第三种组合形式*********/
	//float h1;//u0=-fhan(e1,e2,r,h1);
	//int N1;//跟踪微分器解决速度超调h0=N*h
	/*********第四种组合形式*********/
	//float c;//u0=-fhan(e1,c*e2*e2,r,h1);


}Fhan_Data;


//函数名：ADRC_Init
//作用：初始化
//参数r：微分跟踪器步长，越小跟踪速度越慢，曲线越平缓
//参数h：时间步长
//参数z：线性与非线性边界
//参数p：类似于PID比例环节
//参数d：类似于PID微分环节
void ADRC_Init(Fhan_Data* fhan_Input, float r, float h, float z, float p, float d);
void Fhan_ADRC(Fhan_Data* fhan_Input, float expect_ADRC);
void ADRC_Control(Fhan_Data* fhan_Input, float expect_ADRC, float feedback);

extern Fhan_Data ADRC_dx, ADRC_dy, ADRC_dz;
extern Fhan_Data dz_td;
#endif

