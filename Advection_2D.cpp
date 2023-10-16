#include<stdio.h>
#include<conio.h>
#include<math.h>
#include<string.h>
#define max(x,y) (x>y)?x:y
#define min(x,y) (x<y)?x:y
using namespace std;
int imax=17;
int jmax=17;
double LX=1.0;
double LY=1.0;
double DELX=LX/(imax-2);
double DELY=LY/(jmax-2);
double dt=0.01;
double cp=1.0;
int i,j,k;
int ptime=0;
double rho=1;
//double xx[jmax-1][imax-1],xy[jmax-1][imax-1],x[jmax][imax],y[jmax][imax];
//double U[jmax-1][imax-1],V[jmax-1][imax-1],MX[jmax-1][imax-1],MY[jmax-1][imax-1],HX[jmax-1][imax-1],HY[jmax-1][imax-1],T[jmax][imax],TOLD[jmax][imax],Q[jmax-2][imax-2];
double **xx,**xy,**x,**y,**T,**TOLD,**U,**V,**MX,**MY,**HX,**HY,**Q;
void MAKE_ARRAY();
void SET_GEOMETRY();
void BC();
void IC();
void MASS_FLUX_RATE();
void ENTHALPY();
void UPDATE();
void TEMPERATURE();
void printit(int);
double UNSTEADINESS();
int main(int argc,char *argv[])
{
	double unstead=1;
	int iter=0;
	MAKE_ARRAY();
	SET_GEOMETRY();
	IC();
	BC();
	printit(0);
	while(unstead>0.0001 && iter<=10000)
	{
		iter++;
		unstead=0;
		BC();
		MASS_FLUX_RATE();
		ENTHALPY();
		UPDATE();
		TEMPERATURE();
		if(iter%100 == 0)
		{
			printit(iter);
		}
		unstead=UNSTEADINESS();
		printf("\n Iter=%d Unsteadiness=%f",iter,unstead);
	}
}
void MAKE_ARRAY()
{
	xx=new double*[jmax-1];
	xy=new double*[jmax-1];
	MX=new double*[jmax-1];
	MY=new double*[jmax-1];
	HX=new double*[jmax-1];
	HY=new double*[jmax-1];
	U=new double*[jmax-1];
	V=new double*[jmax-1];
	for(j=0;j<jmax-1;j++)
	{
		xx[j]=new double[imax-1];
		xy[j]=new double[imax-1];
		MX[j]=new double[imax-1];
		MY[j]=new double[imax-1];
		HX[j]=new double[imax-1];
		HY[j]=new double[imax-1];
		U[j]=new double[imax-1];
		V[j]=new double[imax-1];
	}
	x=new double*[jmax];
	y=new double*[jmax];
	T=new double*[jmax];
	TOLD=new double*[jmax];
	for(j=0;j<jmax;j++)
	{
		x[j]=new double[imax];
		y[j]=new double[imax];
		T[j]=new double[imax];
		TOLD[j]=new double[imax];
	}
	Q=new double*[jmax-2];
	for(j=0;j<jmax-2;j++)
	{
		Q[j]=new double[imax-2];
	}
}
void SET_GEOMETRY()
{
	for(j=0;j<jmax-1;j++)
	{
		for(i=0;i<imax-1;i++)
		{
			xx[j][i]=i*DELX;
			xy[j][i]=j*DELY;
		}
	}
	for(j=0;j<jmax;j++)
	{
		for(i=0;i<imax;i++)
		{
			if(i!=0 && i!=1 && i!=imax-1)
			{
				x[j][i]=x[j][i-1]+DELX;
			}
			else if(i==0)
			{
				x[j][i]=0;
			}
			else if(i==1)
			{
				x[j][i]=DELX/2;
			}
			else if(i==imax-1)
			{
				x[j][i]=LX;
			}
		}
	}
	for(j=0;j<jmax;j++)
	{
		for(i=0;i<imax;i++)
		{
			if(j!=0 && j!=1 && j!=jmax-1)
			{
				y[j][i]=y[j-1][i]+DELY;
			}
			else if(j==0)
			{
				y[j][i]=0;
			}
			else if(j==1)
			{
				y[j][i]=DELY/2;
			}
			else if(j==jmax-1)
			{
				y[j][i]=LY;
			}
		}
	}
}
void BC()
{
	for(j=0;j<jmax;j++)
	{
		T[j][0]=100;
		T[j][imax-1]=0;
	}
	for(i=0;i<imax;i++)
	{
		T[0][i]=0;
		T[jmax-1][i]=100;
	}
}
void IC()
{
	for(j=0;j<jmax;j++)
	{
		for(i=0;i<imax;i++)
		{
			T[j][i]=0;
		}
	}
	for(j=0;j<jmax-1;j++)
	{
		for(i=0;i<imax-1;i++)
		{
			U[j][i]=1/sqrt(2);
			V[j][i]=1/sqrt(2);
			MX[j][i]=0;
			MY[j][i]=0;
			HX[j][i]=0;
			HY[j][i]=0;
		}
	}
}
void MASS_FLUX_RATE()
{
	for(j=0;j<jmax-1;j++)
	{
		for(i=0;i<imax-1;i++)
		{
			MX[j][i]=rho*U[j][i];
		}
	}
	for(j=0;j<jmax-1;j++)
	{
		for(i=0;i<imax-1;i++)
		{
			MY[j][i]=rho*V[j][i];
		}
	}
}
void ENTHALPY()
{
	double w1=0.375,w2=0.75,w3=-0.125;
	double mp,mm;
	for(j=1;j<jmax-1;j++)
	{
		for(i=1;i<imax-2;i++)
		{
			mp=max(MX[j][i],0);
			mm=min(MX[j][i],0);
			HX[j][i]=cp*(mp*(w1*T[j][i+1]+w2*T[j][i]+w3*T[j][i-1])+mm*(w1*T[j][i]+w2*T[j][i+1]+w3*T[j][i+2]));
		}
	}
	for(j=1;j<jmax-1;j++)
	{
		mp=max(MX[j][0],0);
		mm=min(MX[j][0],0);
		HX[j][0]=cp*(mp*T[j][0]+mm*T[j][1]);
		mp=max(MX[j][imax-2],0);
		mm=min(MX[j][imax-2],0);
		HX[j][imax-2]=cp*(mp*T[j][imax-2]+mm*T[j][imax-1]);
	}
	for(j=1;j<jmax-2;j++)
	{
		for(i=1;i<imax-1;i++)
		{
			mp=max(MY[j][i],0);
			mm=min(MY[j][i],0);
			HY[j][i]=cp*(mp*(w1*T[j+1][i]+w2*T[j][i]+w3*T[j-1][i])+mm*(w1*T[j][i]+w2*T[j+1][i]+w3*T[j+2][i]));
		}
	}
	for(i=1;i<imax-1;i++)
	{
		mp=max(MY[0][i],0);
		mm=min(MY[0][i],0);
		HY[0][i]=cp*(mp*T[0][i]+mm*T[1][i]);
		mp=max(MY[jmax-2][i],0);
		mm=min(MY[jmax-2][i],0);
		HY[jmax-2][i]=cp*(mp*T[jmax-2][i]+mm*T[jmax-1][i]);
	}
}
void UPDATE()
{
	for(j=0;j<jmax;j++)
	{
		for(i=0;i<imax;i++)
		{
			TOLD[j][i]=T[j][i];
		}
	}
}
void TEMPERATURE()
{
	for(j=1;j<jmax-1;j++)
	{
		for(i=1;i<imax-1;i++)
		{
			Q[j-1][i-1]=((HX[j][i]-HX[j][i-1])*DELY) + ((HY[j][i]-HY[j-1][i])*DELX);
		}
	}
	for(j=1;j<jmax-1;j++)
	{
		for(i=1;i<imax-1;i++)
		{
			T[j][i]=TOLD[j][i]-((dt/(rho*cp*DELX*DELY))*Q[j-1][i-1]);
		}
	}
}
void printit(int iter)
{
	char tem[80]={0};
	int aa;
	char ch[10];
	aa = sprintf(ch, "%d",ptime);
	strcat(tem,ch);
	strcat(tem,".dat");
	FILE *fs;
	fs=fopen(tem,"w");
		for(i=0;i<imax;i++)
		{
			for(j=0;j<jmax;j++)
			{
				if(i==0&&j==0&&k==0)
				{
				fprintf(fs,"VARIABLES = \"X\", \"Y\", \"T\"\n");
				fprintf(fs,"ZONE I=%d, J=%d, F=POINT",imax,jmax);	
				}
				fprintf(fs,"\n%0.6f %0.6f %0.6f ",x[j][i],y[j][i],T[j][i]);
			}
		}
	
	fclose(fs);
	ptime++;
}
double UNSTEADINESS()
{
	double sum=0;
	int cnt=0;
	for(j=1;j<jmax-1;j++)
	{
		for(i=1;i<imax-1;i++)
		{
			cnt++;
			sum+=fabs(T[j][i]-TOLD[j][i]);
		}
	}
	return sqrt(sum/cnt);
}
