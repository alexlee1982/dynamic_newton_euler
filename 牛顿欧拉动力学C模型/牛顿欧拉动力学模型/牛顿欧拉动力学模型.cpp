/*#include<stdio.h>
#include<math.h>
#include "stdlib.h"

 #define pi 3.14159265358979323846

struct m_Vector
{
	double m_shuzu[3];
};

//将数组转换成结构体变量
m_Vector Fuzhi(double shuzu[3])
{
	m_Vector temp;
	temp.m_shuzu[0]=
		shuzu[0];
	temp.m_shuzu[1]=
		shuzu[1];
	temp.m_shuzu[2]=
		shuzu[2];
	return temp;
}

//两向量相乘
m_Vector vectorCross(m_Vector vectorA,m_Vector vectorB)
{
	m_Vector vectorC;

	vectorC.m_shuzu[0]=
		vectorA.m_shuzu[1]*
		vectorB.m_shuzu[2]-
		vectorA.m_shuzu[2]*
		vectorB.m_shuzu[1];
	vectorC.m_shuzu[1]=
		vectorA.m_shuzu[2]*
		vectorB.m_shuzu[0]-
		vectorA.m_shuzu[0]*
		vectorB.m_shuzu[2];
	vectorC.m_shuzu[2]=
		vectorA.m_shuzu[0]*
		vectorB.m_shuzu[1]-
		vectorA.m_shuzu[1]*
		vectorB.m_shuzu[0];

	return vectorC;
}

// 矩阵相乘
m_Vector MultiplyMatrix(double matrix1[3][3],m_Vector vectorA)
{
	m_Vector temp;
	temp.m_shuzu[0]=0;
	temp.m_shuzu[1]=0;
	temp.m_shuzu[2]=0;

	int i,j,k;
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 1; j++)
		{
			for(k = 0; k < 3; k++)
			{
				temp.m_shuzu[i]+= matrix1[i][k]*vectorA.m_shuzu[k];
			}
		}
	}
	return temp;
}

//两向量相加
m_Vector PlusVector(m_Vector vectorA,m_Vector vectorB)
{
	m_Vector temp;
	temp.m_shuzu[0]=0;
	temp.m_shuzu[1]=0;
	temp.m_shuzu[2]=0;

	temp.m_shuzu[0]=
		vectorA.m_shuzu[0]+
		vectorB.m_shuzu[0];
	temp.m_shuzu[1]=
		vectorA.m_shuzu[1]+
		vectorB.m_shuzu[1];
	temp.m_shuzu[2]=
		vectorA.m_shuzu[2]+
		vectorB.m_shuzu[2];

	return temp;
}

//常数与向量相乘
m_Vector ConstPlusVector(double m,m_Vector vectorA)
{
	m_Vector temp;
	temp.m_shuzu[0]=0;
	temp.m_shuzu[1]=0;
	temp.m_shuzu[2]=0;

	temp.m_shuzu[0]=
		m*vectorA.m_shuzu[0];
	temp.m_shuzu[1]=
		m*vectorA.m_shuzu[1];
	temp.m_shuzu[2]=
		m*vectorA.m_shuzu[2];

	return temp;
}

int main()
{
	//////连杆参数单位(米)
	double d1,d3,d4,a1,a2,a3; 
	//double m;
	d1=280e-3;
	d3=0;
	d4=250e-3;
	a1=0;
	a2=d4=250e-3;
	a3=0;
	//m=0.148614701e-3;

	//各连杆质心位置，单位（米）
	int i,j;
	double rc1[3]={0.148614701,0.039458923,-0.057802885};
	m_Vector m_rc1=Fuzhi(rc1);
	double rc2[3]={0.338106843,-0.000269964,0.210197035};
	m_Vector m_rc2=Fuzhi(rc2);
	double rc3[3]={0.082366045,0.053144421,0.028833443};
	m_Vector m_rc3=Fuzhi(rc3);
	double rc4[3]={-0.000000030,0.000000012,-0.365903198};
	m_Vector m_rc4=Fuzhi(rc4);
	double rc5[3]={0.005,0.073505724,0.008};
	m_Vector m_rc5=Fuzhi(rc5);
	double rc6[3]={0.010,0.020,0.280512375};
	m_Vector m_rc6=Fuzhi(rc6);
	//printf("r=%f\n",rc6);

	//%%%%各连杆质心惯性张量，单位kg·m^2（千克*平米）
	double Ic1[3][3]=
	{
		{430428.726553e-6,93198.194264e-6,194134.512662e-6},
		{93198.194264e-6,682005.115471e-6,59214.016610e-6},
		{194134.512662e-6,59214.016610e-6,566666.070984e-6}
	};
	double Ic2[3][3]=
	{
		{67343.901287e-6,-333.945941e-6,28342.482531e-6},
		{-333.945941e-6,1573167.887392e-6,-95.793185e-6},
		{28342.482531e-6,-95.793185e-6,1613465.436738e-6}
	};
	double Ic3[3][3]=
	{
		{189468.304124e-6,63377.667540e-6,-35758.042804e-6},
		{63377.667540e-6,214922.962170e-6,-24220.605180e-6},
		{-35758.042804e-6,-24220.605180e-6,256805.417005e-6}
	};
	double Ic4[3][3]=
	{
		{347447.490456e-6,-0.003695e-6,-0.036357e-6},
		{-0.003695e-6,347413.013567e-6,0.015886e-6},
		{-0.036357e-6,0.015886e-6,7535.874945e-6}
	};
	double Ic5[3][3]=
	{
		25928.997668e-6,500e-6,400e-6,
		500e-6,13460.552429e-6,600e-6,
		400e-6,600e-6,23583.472988e-6
	};
	double Ic6[3][3]=
	{
		11752.881397e-6,-100e-6,-200e-6,
		-100e-6,7994.693551e-6,-300e-6,
		-200e-6,-300e-6,6978.407659e-6
	};

	//各连杆质质量，单位kg（千克）
	double m1=23.107189,m2=11.055208,m3=15.359734,m4=5.188882,m5=3.345372,m6=1.586331;
	//////定义关节空间轨迹
	double	t,c1,dc1,ddc1,c2,dc2,ddc2,c3,dc3,ddc3,c4,dc4,ddc4,c5,dc5,ddc5,c6,dc6,ddc6;

	double m_minTime=0.001;

	FILE *fpt1; 
	fpt1 = fopen("T1.txt","w");//打开文档，写入 
	FILE *fpt2; 
	fpt2 = fopen("T2.txt","w");//打开文档，写入 
	FILE *fpt3; 
	fpt3 = fopen("T3.txt","w");//打开文档，写入 
	FILE *fpt4; 
	fpt4 = fopen("T4.txt","w");//打开文档，写入 
	FILE *fpt5; 
	fpt5 = fopen("T5.txt","w");//打开文档，写入 
	FILE *fpt6; 
	fpt6 = fopen("T6.txt","w");//打开文档，写入 

	for(i=0;i<=5000;i++)
    {
		double t=m_minTime*i;
		double ddc1[3]={0,0,1*pi/180};       //%%定义各关节加速度与时间t的函数
		double ddc2[3]={0,0,1.5*pi/180};
		double ddc3[3]={0,0,2*pi/180};
		double ddc4[3]={0,0,2.5*pi/180};
		double ddc5[3]={0,0,2*pi/180};
		double ddc6[3]={0,0,1*pi/180};

		double dc1[3]={0,0,t*pi/180};         //%%定义各关节速度与时间t的函数
		double dc2[3]={0,0,1.5*t*pi/180};
		double dc3[3]={0,0,2*t*pi/180};
		double dc4[3]={0,0,2.5*t*pi/180};
		double dc5[3]={0,0,2*t*pi/180};
		double dc6[3]={0,0,t*pi/180}; 

		double	 c1=0.5*t*t*pi/180;         //%%定义各关节角度与时间t的函数
		double	 c2=0.75*t*t*pi/180-pi/2;
		double	 c3=t*t*pi/180;
		double	 c4=1.25*t*t*pi/180;
		double	 c5=t*t*pi/180;
		double	 c6=0.5*t*t*pi/180;
		m_Vector m_ddc1=Fuzhi(ddc1);
		m_Vector m_ddc2=Fuzhi(ddc2);
		m_Vector m_ddc3=Fuzhi(ddc3);
		m_Vector m_ddc4=Fuzhi(ddc4);
		m_Vector m_ddc5=Fuzhi(ddc5);
		m_Vector m_ddc6=Fuzhi(ddc6);
		m_Vector m_dc1=Fuzhi(dc1);
		m_Vector m_dc2=Fuzhi(dc2);
		m_Vector m_dc3=Fuzhi(dc3);
		m_Vector m_dc4=Fuzhi(dc4);
		m_Vector m_dc5=Fuzhi(dc5);
		m_Vector m_dc6=Fuzhi(dc6);

		double	 R1[3][3]={cos(c1),sin(c1),0 ,-sin(c1),cos(c1),0,0,0,1};
		double	 R2[3][3]={cos(c2),0,-sin(c2),-sin(c2),0,-cos(c2),0,1,0};
		double	 R3[3][3]={cos(c3),sin(c3),0,-sin(c3),cos(c3),0,0,0,1};
		double	 R4[3][3]={cos(c4),0,-sin(c4),-sin(c4),0,-cos(c4),0,1,0};
		double	 R5[3][3]={cos(c5),0,sin(c5),-sin(c5),0,cos(c5),0,-1,0};
		double	 R6[3][3]={cos(c6),0,-sin(c6),-sin(c6),0,-cos(c6),0,1,0};
		
		double transpose_R1[3][3]={cos(c1),-sin(c1),0,sin(c1),cos(c1),0,0,0,1};
		double transpose_R2[3][3]={cos(c2),-sin(c2),0,0,0,1,-sin(c2),-cos(c2),0};
		double transpose_R3[3][3]={cos(c3),-sin(c3),0,sin(c3),cos(c3),0,0,0,1};
		double transpose_R4[3][3]={cos(c4),-sin(c4),0,0,0,1,-sin(c4),-cos(c4),0};
		double transpose_R5[3][3]={cos(c5),-sin(c5),0,0,0,-1,sin(c5),cos(c5),0};
		double transpose_R6[3][3]={cos(c6),-sin(c6),0,0,0,1,-sin(c6),-cos(c6),0};

		double w0[3]={0,0,0};
		m_Vector m_w0=Fuzhi(w0);
		double dw0[3]={0,0,0};
		m_Vector m_dw0=Fuzhi(dw0);
		double v0[3]={0,0,0};
		m_Vector m_v0=Fuzhi(v0);
		double dv0[3]={0,0,9.80665};
		m_Vector m_dv0=Fuzhi(dv0);
		double f7[3]={0,0,0};
		m_Vector m_f7=Fuzhi(f7);
		double n7[3]={0,0,0};
		m_Vector m_n7=Fuzhi(n7);
		double R7[3][3]={0,0,0,0,0,0,0,0,0};
		double fc7[3]={0,0,0};
		m_Vector m_fc7=Fuzhi(fc7);
		double nc7[3]={0,0,0};
		m_Vector m_nc7=Fuzhi(nc7);
		double r7[3]={0,0,0};
		m_Vector m_r7=Fuzhi(r7);
		double Z[3]={0,0,1};///关节0和关节7参数设置。关节7没有，参数为0，惯性力补偿重力。
		m_Vector m_Z=Fuzhi(Z);
		double r1[3]={0,0,d1};
		m_Vector m_r1=Fuzhi(r1);
		double r2[3]={a1,0,0};
		m_Vector m_r2=Fuzhi(r2);
		double r3[3]={a2,0,d3};
		m_Vector m_r3=Fuzhi(r3);
		double r4[3]={a3,d4,0};
		m_Vector m_r4=Fuzhi(r4);
		double r5[3]={0,0,0};
		m_Vector m_r5=Fuzhi(r5);
		double r6[3]={0,0,0};
		m_Vector m_r6=Fuzhi(r6);

		//第1关节
		m_Vector temp1=MultiplyMatrix(R1,m_w0);
		m_Vector m_w1=PlusVector(temp1,m_dc1);
		//w1=R1*w0+dc1;
		m_Vector m_dw1=PlusVector(PlusVector(MultiplyMatrix(R1,m_dw0),vectorCross(MultiplyMatrix(R1,m_w0),m_dc1)),m_ddc1);
		//dw1=R1*dw0+cross(R1*w0,dc1)+ddc1;
		m_Vector m_v1=MultiplyMatrix(R1,PlusVector(m_v0,vectorCross(m_w0,m_r1)));
		//v1=R1*(v0+cross(w0,r1));
		m_Vector m_dv1=MultiplyMatrix(R1,PlusVector(PlusVector(m_dv0,vectorCross(m_dw0,m_r1)),vectorCross(m_w0,vectorCross(m_w0,m_r1))));
		//dv1=R1*(dv0+cross(dw0,r1)+cross(w0,(cross(w0,r1))));
		m_Vector m_dvc1=PlusVector(PlusVector(vectorCross(m_dw1,m_rc1),vectorCross(m_w1,vectorCross(m_w1,m_rc1))),m_dv1);
		//dvc1=cross(dw1,rc1)+cross(w1,(cross(w1,rc1)))+dv1;
		m_Vector m_fc1=ConstPlusVector(m1,m_dvc1);
		//fc1=m1*dvc1; 
		m_Vector m_nc1=PlusVector(MultiplyMatrix(Ic1,m_dw1),vectorCross(m_w1,MultiplyMatrix(Ic1,m_w1)));
		//nc1=Ic1*dw1+cross(w1,Ic1*w1);

		/////%%第2关节
		m_Vector temp2=MultiplyMatrix(R2,m_w1);
		m_Vector m_w2=PlusVector(temp2,m_dc2);
		//w2=R2*w1+dc2;
		m_Vector m_dw2=PlusVector(PlusVector(MultiplyMatrix(R2,m_dw1),vectorCross(MultiplyMatrix(R2,m_w1),m_dc2)),m_ddc2);
		//dw2=R2*dw1+cross(R2*w1,dc2)+ddc2;
		m_Vector m_v2=MultiplyMatrix(R2,PlusVector(m_v1,vectorCross(m_w1,m_r2)));
		//v2=R2*(v1+cross(w1,r2));
		m_Vector m_dv2=MultiplyMatrix(R2,PlusVector(PlusVector(m_dv1,vectorCross(m_dw1,m_r2)),vectorCross(m_w1,vectorCross(m_w1,m_r2))));
		//dv2=R2*(dv1+cross(dw1,r2)+cross(w1,(cross(w1,r2))));
		m_Vector m_dvc2=PlusVector(PlusVector(vectorCross(m_dw2,m_rc2),vectorCross(m_w2,vectorCross(m_w2,m_rc2))),m_dv2);
		//dvc2=cross(dw2,rc2)+cross(w2,(cross(w2,rc2)))+dv2;
		m_Vector m_fc2=ConstPlusVector(m2,m_dvc2);
		//fc2=m2*dvc2;
		m_Vector m_nc2=PlusVector(MultiplyMatrix(Ic2,m_dw2),vectorCross(m_w2,MultiplyMatrix(Ic2,m_w2)));
		//nc2=Ic2*dw2+cross(w2,Ic2*w2);

		/////%%第3关节
		m_Vector temp3=MultiplyMatrix(R3,m_w2);
		m_Vector m_w3=PlusVector(temp3,m_dc3);
		//w3=R3*w2+dc3;
		m_Vector m_dw3=PlusVector(PlusVector(MultiplyMatrix(R3,m_dw2),vectorCross(MultiplyMatrix(R3,m_w2),m_dc3)),m_ddc3);
		//dw3=R3*dw2+cross(R3*w2,dc3)+ddc3;
		m_Vector m_v3=MultiplyMatrix(R3,PlusVector(m_v2,vectorCross(m_w2,m_r3)));
		//v3=R3*(v2+cross(w2,r3));
		m_Vector m_dv3=MultiplyMatrix(R3,PlusVector(PlusVector(m_dv2,vectorCross(m_dw2,m_r3)),vectorCross(m_w2,vectorCross(m_w2,m_r3))));
		//dv3=R3*(dv2+cross(dw2,r3)+cross(w2,(cross(w2,r3))));
		m_Vector m_dvc3=PlusVector(PlusVector(vectorCross(m_dw3,m_rc3),vectorCross(m_w3,vectorCross(m_w3,m_rc3))),m_dv3);
		//dvc3=cross(dw3,rc3)+cross(w3,(cross(w3,rc3)))+dv3;
		m_Vector m_fc3=ConstPlusVector(m3,m_dvc3);
		//fc3=m3*dvc3;
		m_Vector m_nc3=PlusVector(MultiplyMatrix(Ic3,m_dw3),vectorCross(m_w3,MultiplyMatrix(Ic3,m_w3)));
		//nc3=Ic3*dw3+cross(w3,Ic3*w3);

		/////%%第4关节
		m_Vector temp4=MultiplyMatrix(R4,m_w3);
		m_Vector m_w4=PlusVector(temp4,m_dc4);
		m_Vector m_dw4=PlusVector(PlusVector(MultiplyMatrix(R4,m_dw3),vectorCross(MultiplyMatrix(R4,m_w3),m_dc4)),m_ddc4);
		m_Vector m_v4=MultiplyMatrix(R4,PlusVector(m_v3,vectorCross(m_w3,m_r4)));
		m_Vector m_dv4=MultiplyMatrix(R4,PlusVector(PlusVector(m_dv3,vectorCross(m_dw3,m_r4)),vectorCross(m_w3,vectorCross(m_w3,m_r4))));
		m_Vector m_dvc4=PlusVector(PlusVector(vectorCross(m_dw4,m_rc4),vectorCross(m_w4,vectorCross(m_w4,m_rc4))),m_dv4);
		m_Vector m_fc4=ConstPlusVector(m4,m_dvc4);
		m_Vector m_nc4=PlusVector(MultiplyMatrix(Ic4,m_dw4),vectorCross(m_w4,MultiplyMatrix(Ic4,m_w4)));

		/////%%第5关节
		m_Vector temp5=MultiplyMatrix(R5,m_w4);
		m_Vector m_w5=PlusVector(temp5,m_dc5);
		m_Vector m_dw5=PlusVector(PlusVector(MultiplyMatrix(R5,m_dw4),vectorCross(MultiplyMatrix(R5,m_w4),m_dc5)),m_ddc5);
		m_Vector m_v5=MultiplyMatrix(R5,PlusVector(m_v4,vectorCross(m_w4,m_r5)));
		m_Vector m_dv5=MultiplyMatrix(R5,PlusVector(PlusVector(m_dv4,vectorCross(m_dw4,m_r5)),vectorCross(m_w4,vectorCross(m_w4,m_r5))));
		m_Vector m_dvc5=PlusVector(PlusVector(vectorCross(m_dw5,m_rc5),vectorCross(m_w5,vectorCross(m_w5,m_rc5))),m_dv5);
		m_Vector m_fc5=ConstPlusVector(m5,m_dvc5);
		m_Vector m_nc5=PlusVector(MultiplyMatrix(Ic5,m_dw5),vectorCross(m_w5,MultiplyMatrix(Ic5,m_w5)));

		/////%%第6关节
		m_Vector temp6=MultiplyMatrix(R6,m_w5);
		m_Vector m_w6=PlusVector(temp6,m_dc6);
		m_Vector m_dw6=PlusVector(PlusVector(MultiplyMatrix(R6,m_dw5),vectorCross(MultiplyMatrix(R6,m_w5),m_dc6)),m_ddc6);
		m_Vector m_v6=MultiplyMatrix(R6,PlusVector(m_v5,vectorCross(m_w5,m_r6)));
		m_Vector m_dv6=MultiplyMatrix(R6,PlusVector(PlusVector(m_dv5,vectorCross(m_dw5,m_r6)),vectorCross(m_w5,vectorCross(m_w5,m_r6))));
		m_Vector m_dvc6=PlusVector(PlusVector(vectorCross(m_dw6,m_rc6),vectorCross(m_w6,vectorCross(m_w6,m_rc6))),m_dv6);
		m_Vector m_fc6=ConstPlusVector(m6,m_dvc6);
		m_Vector m_nc6=PlusVector(MultiplyMatrix(Ic6,m_dw6),vectorCross(m_w6,MultiplyMatrix(Ic6,m_w6)));
		
		/////内推：i：6→1.
		////%%%第6关节
		m_Vector m_f6=PlusVector(MultiplyMatrix(R7,m_f7),m_fc6);
		//f6=R7*f7+fc6;
		m_Vector m_n6=PlusVector(PlusVector(PlusVector(m_nc6,MultiplyMatrix(R7,m_n7)),vectorCross(m_rc6,m_fc6)),vectorCross(m_r7,MultiplyMatrix(R7,m_f7)));
		//n6=nc6+R7*n7+cross(rc6,fc6)+cross(r7,R7*f7);
		double t6=m_n6.m_shuzu[2]*1000;

		////%%%第5关节
		m_Vector m_f5=PlusVector(MultiplyMatrix(transpose_R6,m_f6),m_fc5);
		m_Vector m_n5=PlusVector(PlusVector(PlusVector(m_nc5,MultiplyMatrix(transpose_R6,m_n6)),vectorCross(m_rc5,m_fc5)),vectorCross(m_r6,MultiplyMatrix(transpose_R6,m_f6)));
		//n5=nc5+R6'*n6+cross(rc5,fc5)+cross(r6,R6'*f6);
		double t5=m_n5.m_shuzu[2]*1000;
	
		////%%%第4关节
		m_Vector m_f4=PlusVector(MultiplyMatrix(transpose_R5,m_f5),m_fc4);
		m_Vector m_n4=PlusVector(PlusVector(PlusVector(m_nc4,MultiplyMatrix(transpose_R5,m_n5)),vectorCross(m_rc4,m_fc4)),vectorCross(m_r5,MultiplyMatrix(transpose_R5,m_f5)));
		double t4=m_n4.m_shuzu[2]*1000;

		////%%%第3关节
		m_Vector m_f3=PlusVector(MultiplyMatrix(transpose_R4,m_f4),m_fc3);
		m_Vector m_n3=PlusVector(PlusVector(PlusVector(m_nc3,MultiplyMatrix(transpose_R4,m_n4)),vectorCross(m_rc3,m_fc3)),vectorCross(m_r4,MultiplyMatrix(transpose_R4,m_f4)));
		double t3=m_n3.m_shuzu[2]*1000;

		////%%%第2关节
		m_Vector m_f2=PlusVector(MultiplyMatrix(transpose_R3,m_f3),m_fc2);
		m_Vector m_n2=PlusVector(PlusVector(PlusVector(m_nc2,MultiplyMatrix(transpose_R3,m_n3)),vectorCross(m_rc2,m_fc2)),vectorCross(m_r3,MultiplyMatrix(transpose_R3,m_f3)));
		double t2=m_n2.m_shuzu[2]*1000;

		////%%%第1关节
		m_Vector m_f1=PlusVector(MultiplyMatrix(transpose_R2,m_f2),m_fc1);
		m_Vector m_n1=PlusVector(PlusVector(PlusVector(m_nc1,MultiplyMatrix(transpose_R2,m_n2)),vectorCross(m_rc1,m_fc1)),vectorCross(m_r2,MultiplyMatrix(transpose_R2,m_f2)));
		double t1=m_n1.m_shuzu[2]*1000;
		fprintf(fpt1,"%8.8f\r\n",t1); 
		fprintf(fpt2,"%8.8f\r\n",t2); 
		fprintf(fpt3,"%8.8f\r\n",t3); 
		fprintf(fpt4,"%8.8f\r\n",t4); 
		fprintf(fpt5,"%8.8f\r\n",t5); 
		fprintf(fpt6,"%8.8f\r\n",t6); 
	}
	fclose(fpt1); 
	fclose(fpt2);
	fclose(fpt3); 
	fclose(fpt4); 
	fclose(fpt5); 
	fclose(fpt6); 
}*/
#include<stdio.h>
#include<math.h>
#include "stdlib.h"

 #define pi 3.14159265358979323846264338327950

struct m_Vector
{
	double m_shuzu[3];
};

//将数组转换成结构体变量
m_Vector Fuzhi(double shuzu[3])
{
	m_Vector temp;
	temp.m_shuzu[0]=
		shuzu[0];
	temp.m_shuzu[1]=
		shuzu[1];
	temp.m_shuzu[2]=
		shuzu[2];
	return temp;
}

//两向量相乘
m_Vector vectorCross(m_Vector vectorA,m_Vector vectorB)
{
	m_Vector vectorC;

	vectorC.m_shuzu[0]=
		vectorA.m_shuzu[1]*
		vectorB.m_shuzu[2]-
		vectorA.m_shuzu[2]*
		vectorB.m_shuzu[1];
	vectorC.m_shuzu[1]=
		vectorA.m_shuzu[2]*
		vectorB.m_shuzu[0]-
		vectorA.m_shuzu[0]*
		vectorB.m_shuzu[2];
	vectorC.m_shuzu[2]=
		vectorA.m_shuzu[0]*
		vectorB.m_shuzu[1]-
		vectorA.m_shuzu[1]*
		vectorB.m_shuzu[0];

	return vectorC;
}

// 矩阵相乘
m_Vector MultiplyMatrix(double matrix1[3][3],m_Vector vectorA)
{
	m_Vector temp;
	temp.m_shuzu[0]=0;
	temp.m_shuzu[1]=0;
	temp.m_shuzu[2]=0;

	int i,j,k;
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 1; j++)
		{
			for(k = 0; k < 3; k++)
			{
				temp.m_shuzu[i]+= matrix1[i][k]*vectorA.m_shuzu[k];
			}
		}
	}
	return temp;
}

//两向量相加
m_Vector PlusVector(m_Vector vectorA,m_Vector vectorB)
{
	m_Vector temp;
	temp.m_shuzu[0]=0;
	temp.m_shuzu[1]=0;
	temp.m_shuzu[2]=0;

	temp.m_shuzu[0]=
		vectorA.m_shuzu[0]+
		vectorB.m_shuzu[0];
	temp.m_shuzu[1]=
		vectorA.m_shuzu[1]+
		vectorB.m_shuzu[1];
	temp.m_shuzu[2]=
		vectorA.m_shuzu[2]+
		vectorB.m_shuzu[2];

	return temp;
}

//常数与向量相乘
m_Vector ConstPlusVector(double m,m_Vector vectorA)
{
	m_Vector temp;
	temp.m_shuzu[0]=0;
	temp.m_shuzu[1]=0;
	temp.m_shuzu[2]=0;

	temp.m_shuzu[0]=
		m*vectorA.m_shuzu[0];
	temp.m_shuzu[1]=
		m*vectorA.m_shuzu[1];
	temp.m_shuzu[2]=
		m*vectorA.m_shuzu[2];

	return temp;
}

int main()
{
	//////连杆参数单位(米)
	double d1,d3,d4,a1,a2,a3; 
	//double m;
	d1=280e-3;
	d3=0;
	d4=1013.5e-3;
	a1=220e-3;
	a2=900e-3;
	a3=160e-3;
	//m=0.148614701e-3;

	//各连杆质心位置，单位（米）
	int i,j;
	double rc1[3]={ 0.01944729909, 0.02980477858, -0.22180483830 };
	m_Vector m_rc1=Fuzhi(rc1);
	double rc2[3]={ 0.31860220304, -0.00, 0.16280837626 };
	m_Vector m_rc2=Fuzhi(rc2);
	double rc3[3]={ 0.06390342865, 0.07151117703, 0.08934224425 };
	m_Vector m_rc3=Fuzhi(rc3);
	double rc4[3]={ 0, 0, -0.39693069721 };
	m_Vector m_rc4=Fuzhi(rc4);
	double rc5[3]={ 0, 0.00739728381, 0 };
	m_Vector m_rc5=Fuzhi(rc5);
	double rc6[3]={ 0, 0, 0.12481798715 };
	m_Vector m_rc6=Fuzhi(rc6);
	//printf("r=%f\n",rc6);

	//%%%%各连杆质心惯性张量，单位kg·m^2（千克*平米）
	double Ic1[3][3]=
{ 1.16097741304095, 0.05138023956666, 0.23052213875993,
			                 0.05138023956666, 1.44440474511337, 0.21098308091826,
			                 0.23052213875993, 0.21098308091826, 0.94785365138672 };
	double Ic2[3][3]=
{ 0.53560927063693, -0.00, 0.29794836700128,
			                 -0.00, 6.87733401521884, -0.00,
			                 0.29794836700128, -0.00, 7.12086740624843 };

double Ic3[3][3] = { 0.27814643940732, 0.09002959963086, -0.07998135854618,
			                 0.09002959963086, 0.31382570555938, -0.10285183833350,
			                 -0.07998135854618, -0.10285183833350, 0.31811802288226 };
		double Ic4[3][3] = { 1.20476946375220, 0.0, 0.0,
			                 0.0, 1.20386579394007, 0.0,
			                 0.0, 0.0, 0.04653310536757 };
		double Ic5[3][3] = { 0.00498517957371, 0.00000000001373, 0.00000000001217,
			                 0.00000000001373, 0.00422321290017, 0.00000000,
			                 0.00000000001217, 0.00000000, 0.00269095479625 };
		double Ic6[3][3] = { 0.00478590437385, 0.000, 0.000,
			                 0.000, 0.00478590437385, 0.000,
			                 0.000, 0.000, 0.00146759697162 };
	
	//各连杆质质量，单位kg（千克）
	double m1=58.45804846,m2=76.70623808,m3=18.10297679,m4=23.87605416,m5=2.17295423,m6=1.98061709;
	//////定义关节空间轨迹
	double	t,c1,dc1,ddc1,c2,dc2,ddc2,c3,dc3,ddc3,c4,dc4,ddc4,c5,dc5,ddc5,c6,dc6,ddc6;

	double m_minTime=0.5;

	FILE *fpt1; 
	fpt1 = fopen("T1.txt","w");//打开文档，写入 
	FILE *fpt2; 
	fpt2 = fopen("T2.txt","w");//打开文档，写入 
	FILE *fpt3; 
	fpt3 = fopen("T3.txt","w");//打开文档，写入 
	FILE *fpt4; 
	fpt4 = fopen("T4.txt","w");//打开文档，写入 
	FILE *fpt5; 
	fpt5 = fopen("T5.txt","w");//打开文档，写入 
	FILE *fpt6; 
	fpt6 = fopen("T6.txt","w");//打开文档，写入 

	for(i=0;i<=10;i++)
    {
		double t=m_minTime*i;
		double ddc1[3]={0,0,1*pi/180};       //%%定义各关节加速度与时间t的函数
		double ddc2[3]={0,0,1.5*pi/180};
		double ddc3[3]={0,0,2*pi/180};
		double ddc4[3]={0,0,2*pi/180};
		double ddc5[3]={0,0,-1.5*pi/180};
		double ddc6[3]={0,0,1*pi/180};

		double dc1[3]={0,0,t*pi/180};         //%%定义各关节速度与时间t的函数
		double dc2[3]={0,0,1.5*t*pi/180};
		double dc3[3]={0,0,2*t*pi/180};
		double dc4[3]={0,0,2*t*pi/180};
		double dc5[3]={0,0,-1.5*t*pi/180};
		double dc6[3]={0,0,t*pi/180}; 

		double	 c1=0.5*t*t*pi/180;         //%%定义各关节角度与时间t的函数
		double	 c2=0.75*t*t*pi/180-pi/2;
		double	 c3=t*t*pi/180;
		double	 c4=1*t*t*pi/180;
		double	 c5=-0.75*t*t*pi/180+pi/4;
		double	 c6=0.5*t*t*pi/180;
		m_Vector m_ddc1=Fuzhi(ddc1);
		m_Vector m_ddc2=Fuzhi(ddc2);
		m_Vector m_ddc3=Fuzhi(ddc3);
		m_Vector m_ddc4=Fuzhi(ddc4);
		m_Vector m_ddc5=Fuzhi(ddc5);
		m_Vector m_ddc6=Fuzhi(ddc6);
		m_Vector m_dc1=Fuzhi(dc1);
		m_Vector m_dc2=Fuzhi(dc2);
		m_Vector m_dc3=Fuzhi(dc3);
		m_Vector m_dc4=Fuzhi(dc4);
		m_Vector m_dc5=Fuzhi(dc5);
		m_Vector m_dc6=Fuzhi(dc6);

		double	 R1[3][3]={cos(c1),sin(c1),0 ,-sin(c1),cos(c1),0,0,0,1};
		double	 R2[3][3]={cos(c2),0,-sin(c2),-sin(c2),0,-cos(c2),0,1,0};
		double	 R3[3][3]={cos(c3),sin(c3),0,-sin(c3),cos(c3),0,0,0,1};
		double	 R4[3][3]={cos(c4),0,-sin(c4),-sin(c4),0,-cos(c4),0,1,0};
		double	 R5[3][3]={cos(c5),0,sin(c5),-sin(c5),0,cos(c5),0,-1,0};
		double	 R6[3][3]={cos(c6),0,-sin(c6),-sin(c6),0,-cos(c6),0,1,0};
		
		double transpose_R1[3][3]={cos(c1),-sin(c1),0,sin(c1),cos(c1),0,0,0,1};
		double transpose_R2[3][3]={cos(c2),-sin(c2),0,0,0,1,-sin(c2),-cos(c2),0};
		double transpose_R3[3][3]={cos(c3),-sin(c3),0,sin(c3),cos(c3),0,0,0,1};
		double transpose_R4[3][3]={cos(c4),-sin(c4),0,0,0,1,-sin(c4),-cos(c4),0};
		double transpose_R5[3][3]={cos(c5),-sin(c5),0,0,0,-1,sin(c5),cos(c5),0};
		double transpose_R6[3][3]={cos(c6),-sin(c6),0,0,0,1,-sin(c6),-cos(c6),0};

		double w0[3]={0,0,0};
		m_Vector m_w0=Fuzhi(w0);
		double dw0[3]={0,0,0};
		m_Vector m_dw0=Fuzhi(dw0);
		double v0[3]={0,0,0};
		m_Vector m_v0=Fuzhi(v0);
		double dv0[3]={0,0,9.80665};
		m_Vector m_dv0=Fuzhi(dv0);
		double f7[3]={0,0,0};
		m_Vector m_f7=Fuzhi(f7);
		double n7[3]={0,0,0};
		m_Vector m_n7=Fuzhi(n7);
		double R7[3][3]={0,0,0,0,0,0,0,0,0};
		double fc7[3]={0,0,0};
		m_Vector m_fc7=Fuzhi(fc7);
		double nc7[3]={0,0,0};
		m_Vector m_nc7=Fuzhi(nc7);
		double r7[3]={0,0,0};
		m_Vector m_r7=Fuzhi(r7);
		double Z[3]={0,0,1};///关节0和关节7参数设置。关节7没有，参数为0，惯性力补偿重力。
		m_Vector m_Z=Fuzhi(Z);
		double r1[3]={0,0,d1};
		m_Vector m_r1=Fuzhi(r1);
		double r2[3]={a1,0,0};
		m_Vector m_r2=Fuzhi(r2);
		double r3[3]={a2,0,d3};
		m_Vector m_r3=Fuzhi(r3);
		double r4[3]={a3,d4,0};
		m_Vector m_r4=Fuzhi(r4);
		double r5[3]={0,0,0};
		m_Vector m_r5=Fuzhi(r5);
		double r6[3]={0,0,0};
		m_Vector m_r6=Fuzhi(r6);

		//第1关节
		m_Vector temp1=MultiplyMatrix(R1,m_w0);
		m_Vector m_w1=PlusVector(temp1,m_dc1);
		//w1=R1*w0+dc1;
		m_Vector m_dw1=PlusVector(PlusVector(MultiplyMatrix(R1,m_dw0),vectorCross(MultiplyMatrix(R1,m_w0),m_dc1)),m_ddc1);
		//dw1=R1*dw0+cross(R1*w0,dc1)+ddc1;
		m_Vector m_v1=MultiplyMatrix(R1,PlusVector(m_v0,vectorCross(m_w0,m_r1)));
		//v1=R1*(v0+cross(w0,r1));
		m_Vector m_dv1=MultiplyMatrix(R1,PlusVector(PlusVector(m_dv0,vectorCross(m_dw0,m_r1)),vectorCross(m_w0,vectorCross(m_w0,m_r1))));
		//dv1=R1*(dv0+cross(dw0,r1)+cross(w0,(cross(w0,r1))));
		m_Vector m_dvc1=PlusVector(PlusVector(vectorCross(m_dw1,m_rc1),vectorCross(m_w1,vectorCross(m_w1,m_rc1))),m_dv1);
		//dvc1=cross(dw1,rc1)+cross(w1,(cross(w1,rc1)))+dv1;
		m_Vector m_fc1=ConstPlusVector(m1,m_dvc1);
		//fc1=m1*dvc1; 
		m_Vector m_nc1=PlusVector(MultiplyMatrix(Ic1,m_dw1),vectorCross(m_w1,MultiplyMatrix(Ic1,m_w1)));
		//nc1=Ic1*dw1+cross(w1,Ic1*w1);

		/////%%第2关节
		m_Vector temp2=MultiplyMatrix(R2,m_w1);
		m_Vector m_w2=PlusVector(temp2,m_dc2);
		//w2=R2*w1+dc2;
		m_Vector m_dw2=PlusVector(PlusVector(MultiplyMatrix(R2,m_dw1),vectorCross(MultiplyMatrix(R2,m_w1),m_dc2)),m_ddc2);
		//dw2=R2*dw1+cross(R2*w1,dc2)+ddc2;
		m_Vector m_v2=MultiplyMatrix(R2,PlusVector(m_v1,vectorCross(m_w1,m_r2)));
		//v2=R2*(v1+cross(w1,r2));
		m_Vector m_dv2=MultiplyMatrix(R2,PlusVector(PlusVector(m_dv1,vectorCross(m_dw1,m_r2)),vectorCross(m_w1,vectorCross(m_w1,m_r2))));
		//dv2=R2*(dv1+cross(dw1,r2)+cross(w1,(cross(w1,r2))));
		m_Vector m_dvc2=PlusVector(PlusVector(vectorCross(m_dw2,m_rc2),vectorCross(m_w2,vectorCross(m_w2,m_rc2))),m_dv2);
		//dvc2=cross(dw2,rc2)+cross(w2,(cross(w2,rc2)))+dv2;
		m_Vector m_fc2=ConstPlusVector(m2,m_dvc2);
		//fc2=m2*dvc2;
		m_Vector m_nc2=PlusVector(MultiplyMatrix(Ic2,m_dw2),vectorCross(m_w2,MultiplyMatrix(Ic2,m_w2)));
		//nc2=Ic2*dw2+cross(w2,Ic2*w2);

		/////%%第3关节
		m_Vector temp3=MultiplyMatrix(R3,m_w2);
		m_Vector m_w3=PlusVector(temp3,m_dc3);
		//w3=R3*w2+dc3;
		m_Vector m_dw3=PlusVector(PlusVector(MultiplyMatrix(R3,m_dw2),vectorCross(MultiplyMatrix(R3,m_w2),m_dc3)),m_ddc3);
		//dw3=R3*dw2+cross(R3*w2,dc3)+ddc3;
		m_Vector m_v3=MultiplyMatrix(R3,PlusVector(m_v2,vectorCross(m_w2,m_r3)));
		//v3=R3*(v2+cross(w2,r3));
		m_Vector m_dv3=MultiplyMatrix(R3,PlusVector(PlusVector(m_dv2,vectorCross(m_dw2,m_r3)),vectorCross(m_w2,vectorCross(m_w2,m_r3))));
		//dv3=R3*(dv2+cross(dw2,r3)+cross(w2,(cross(w2,r3))));
		m_Vector m_dvc3=PlusVector(PlusVector(vectorCross(m_dw3,m_rc3),vectorCross(m_w3,vectorCross(m_w3,m_rc3))),m_dv3);
		//dvc3=cross(dw3,rc3)+cross(w3,(cross(w3,rc3)))+dv3;
		m_Vector m_fc3=ConstPlusVector(m3,m_dvc3);
		//fc3=m3*dvc3;
		m_Vector m_nc3=PlusVector(MultiplyMatrix(Ic3,m_dw3),vectorCross(m_w3,MultiplyMatrix(Ic3,m_w3)));
		//nc3=Ic3*dw3+cross(w3,Ic3*w3);

		/////%%第4关节
		m_Vector temp4=MultiplyMatrix(R4,m_w3);
		m_Vector m_w4=PlusVector(temp4,m_dc4);
		m_Vector m_dw4=PlusVector(PlusVector(MultiplyMatrix(R4,m_dw3),vectorCross(MultiplyMatrix(R4,m_w3),m_dc4)),m_ddc4);
		m_Vector m_v4=MultiplyMatrix(R4,PlusVector(m_v3,vectorCross(m_w3,m_r4)));
		m_Vector m_dv4=MultiplyMatrix(R4,PlusVector(PlusVector(m_dv3,vectorCross(m_dw3,m_r4)),vectorCross(m_w3,vectorCross(m_w3,m_r4))));
		m_Vector m_dvc4=PlusVector(PlusVector(vectorCross(m_dw4,m_rc4),vectorCross(m_w4,vectorCross(m_w4,m_rc4))),m_dv4);
		m_Vector m_fc4=ConstPlusVector(m4,m_dvc4);
		m_Vector m_nc4=PlusVector(MultiplyMatrix(Ic4,m_dw4),vectorCross(m_w4,MultiplyMatrix(Ic4,m_w4)));

		/////%%第5关节
		m_Vector temp5=MultiplyMatrix(R5,m_w4);
		m_Vector m_w5=PlusVector(temp5,m_dc5);
		m_Vector m_dw5=PlusVector(PlusVector(MultiplyMatrix(R5,m_dw4),vectorCross(MultiplyMatrix(R5,m_w4),m_dc5)),m_ddc5);
		m_Vector m_v5=MultiplyMatrix(R5,PlusVector(m_v4,vectorCross(m_w4,m_r5)));
		m_Vector m_dv5=MultiplyMatrix(R5,PlusVector(PlusVector(m_dv4,vectorCross(m_dw4,m_r5)),vectorCross(m_w4,vectorCross(m_w4,m_r5))));
		m_Vector m_dvc5=PlusVector(PlusVector(vectorCross(m_dw5,m_rc5),vectorCross(m_w5,vectorCross(m_w5,m_rc5))),m_dv5);
		m_Vector m_fc5=ConstPlusVector(m5,m_dvc5);
		m_Vector m_nc5=PlusVector(MultiplyMatrix(Ic5,m_dw5),vectorCross(m_w5,MultiplyMatrix(Ic5,m_w5)));

		/////%%第6关节
		m_Vector temp6=MultiplyMatrix(R6,m_w5);
		m_Vector m_w6=PlusVector(temp6,m_dc6);
		m_Vector m_dw6=PlusVector(PlusVector(MultiplyMatrix(R6,m_dw5),vectorCross(MultiplyMatrix(R6,m_w5),m_dc6)),m_ddc6);
		m_Vector m_v6=MultiplyMatrix(R6,PlusVector(m_v5,vectorCross(m_w5,m_r6)));
		m_Vector m_dv6=MultiplyMatrix(R6,PlusVector(PlusVector(m_dv5,vectorCross(m_dw5,m_r6)),vectorCross(m_w5,vectorCross(m_w5,m_r6))));
		m_Vector m_dvc6=PlusVector(PlusVector(vectorCross(m_dw6,m_rc6),vectorCross(m_w6,vectorCross(m_w6,m_rc6))),m_dv6);
		m_Vector m_fc6=ConstPlusVector(m6,m_dvc6);
		m_Vector m_nc6=PlusVector(MultiplyMatrix(Ic6,m_dw6),vectorCross(m_w6,MultiplyMatrix(Ic6,m_w6)));
		
		/////内推：i：6→1.
		////%%%第6关节
		m_Vector m_f6=PlusVector(MultiplyMatrix(R7,m_f7),m_fc6);
		//f6=R7*f7+fc6;
		m_Vector m_n6=PlusVector(PlusVector(PlusVector(m_nc6,MultiplyMatrix(R7,m_n7)),vectorCross(m_rc6,m_fc6)),vectorCross(m_r7,MultiplyMatrix(R7,m_f7)));
		//n6=nc6+R7*n7+cross(rc6,fc6)+cross(r7,R7*f7);
		double t6=m_n6.m_shuzu[2]*1000;

		////%%%第5关节
		m_Vector m_f5=PlusVector(MultiplyMatrix(transpose_R6,m_f6),m_fc5);
		m_Vector m_n5=PlusVector(PlusVector(PlusVector(m_nc5,MultiplyMatrix(transpose_R6,m_n6)),vectorCross(m_rc5,m_fc5)),vectorCross(m_r6,MultiplyMatrix(transpose_R6,m_f6)));
		//n5=nc5+R6'*n6+cross(rc5,fc5)+cross(r6,R6'*f6);
		double t5=m_n5.m_shuzu[2]*1000;
	
		////%%%第4关节
		m_Vector m_f4=PlusVector(MultiplyMatrix(transpose_R5,m_f5),m_fc4);
		m_Vector m_n4=PlusVector(PlusVector(PlusVector(m_nc4,MultiplyMatrix(transpose_R5,m_n5)),vectorCross(m_rc4,m_fc4)),vectorCross(m_r5,MultiplyMatrix(transpose_R5,m_f5)));
		double t4=m_n4.m_shuzu[2]*1000;

		////%%%第3关节
		m_Vector m_f3=PlusVector(MultiplyMatrix(transpose_R4,m_f4),m_fc3);
		m_Vector m_n3=PlusVector(PlusVector(PlusVector(m_nc3,MultiplyMatrix(transpose_R4,m_n4)),vectorCross(m_rc3,m_fc3)),vectorCross(m_r4,MultiplyMatrix(transpose_R4,m_f4)));
		double t3=m_n3.m_shuzu[2]*1000;

		////%%%第2关节
		m_Vector m_f2=PlusVector(MultiplyMatrix(transpose_R3,m_f3),m_fc2);
		m_Vector m_n2=PlusVector(PlusVector(PlusVector(m_nc2,MultiplyMatrix(transpose_R3,m_n3)),vectorCross(m_rc2,m_fc2)),vectorCross(m_r3,MultiplyMatrix(transpose_R3,m_f3)));
		double t2=m_n2.m_shuzu[2]*1000;

		////%%%第1关节
		m_Vector m_f1=PlusVector(MultiplyMatrix(transpose_R2,m_f2),m_fc1);
		m_Vector m_n1=PlusVector(PlusVector(PlusVector(m_nc1,MultiplyMatrix(transpose_R2,m_n2)),vectorCross(m_rc1,m_fc1)),vectorCross(m_r2,MultiplyMatrix(transpose_R2,m_f2)));
		double t1=m_n1.m_shuzu[2]*1000;
		fprintf(fpt1,"%8.15f\r\n",t1); 
		fprintf(fpt2,"%8.15f\r\n",t2); 
		fprintf(fpt3,"%8.15f\r\n",t3); 
		fprintf(fpt4,"%8.15f\r\n",t4); 
		fprintf(fpt5,"%8.15f\r\n",t5); 
		fprintf(fpt6,"%8.15f\r\n",t6); 
	}
	fclose(fpt1); 
	fclose(fpt2);
	fclose(fpt3); 
	fclose(fpt4); 
	fclose(fpt5); 
	fclose(fpt6); 
}


