#include"ADRC.h"
#include<iostream>
using namespace std;
void main()
{
	Fhan_Data x;
	float y = 0;
	ADRC_Init(&x, 1800, 0.0025, 5, 5, 0.01);
	ADRC_Control(&x, 0, y);
	for (int i = 0; i < 100; i++)
	{
		ADRC_Control(&x,30,y);
		y += 0.1f * x.u;
		cout  << x.e1 << ' ' <<  x.e2  << ' ' << y << ' ' << x.u0 << endl;
	}
}