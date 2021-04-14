#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define N 256

//declare functions
void get_surface(double [],double [],double [],double []);
void gauss(double[], double[], double[], int);
double arr_sum(double[],int);
double random(int ,int);

int main()
{
    double xs[5] = {-10,10,-10,10,0};
    double ys[5] = {10,10,-10,-10,0};
    double zs[5] = {3,3,3,3,0};
    double surface[N*N] = {0};
	
    get_surface(xs,ys,zs,surface);
	
    for(int i=0;i<256;i++)
    {
    	for(int j=0;j<256;j++)
    	{	
    		printf("%f  ",surface[i*256+j]);
		}
		printf("\n\n");
	}
		
	return 0;
}

void get_surface(double xs[],double ys[],double zs[],double surface[])
{
       
    int n = 5;   
    double tmp[n] = {0};

    double sigma_x1 = arr_sum(xs,n);
    double sigma_y1 = arr_sum(ys,n);
    double sigma_z1 = arr_sum(zs,5); 

    //sigma_x2
    for(int i=0;i<n;i++)	
	{
		tmp[i] = pow(xs[i],2);
	}
	double sigma_x2 = arr_sum(tmp,n);
	
	
	//sigma_y2
	for(int i=0;i<n;i++)	
	{
		tmp[i] = pow(ys[i],2);
	}
	double sigma_y2 = arr_sum(tmp,n);
	

   	//sigma_x3 
	for(int i=0;i<n;i++)	
	{
		tmp[i] = pow(xs[i],3);
	}
	double sigma_x3 = arr_sum(tmp,n);


    	//sigma_y3 
	for(int i=0;i<n;i++)	
	{
		tmp[i] = pow(ys[i],3);
	}
	double sigma_y3 = arr_sum(tmp,n);	


    //sigma_x4
    for(int i=0;i<n;i++)
    {
    	tmp[i] = pow(xs[i],4);
	}
	double sigma_x4 = arr_sum(tmp,n);
	
	
    //sigma_y4
    for(int i=0;i<n;i++)
    {
    	tmp[i] = pow(ys[i],4);
	}
	double sigma_y4 = arr_sum(tmp,n);


    //sigma_x1y1
    for(int i=0;i<n;i++)
    {
    	tmp[i] = xs[i]*ys[i];
	}
    double sigma_x1y1 = arr_sum(tmp,n);
    
    
    //sigma_x1y2
    for(int i=0;i<n;i++)
    {
    	tmp[i] = xs[i]*ys[i]*ys[i];
	}
    double sigma_x1y2 = arr_sum(tmp,n);
    
    
    //sigma_x1y3
    for(int i=0;i<n;i++)
    {
    	tmp[i] = xs[i]*ys[i]*ys[i]*ys[i];
	}
    double sigma_x1y3 = arr_sum(tmp,n);
    
    
    //sigma_x2y1
    for(int i=0;i<n;i++)
    {
    	tmp[i] = xs[i]*xs[i]*ys[i];
	}
    double sigma_x2y1 = arr_sum(tmp,n);    
    
    
    //sigma_x3y1
    for(int i=0;i<n;i++)
    {
    	tmp[i] = xs[i]*xs[i]*xs[i]*ys[i];
	}
    double sigma_x3y1 = arr_sum(tmp,n);     
    
    
    //sigma_x2y2
    for(int i=0;i<n;i++)
    {
    	tmp[i] = xs[i]*xs[i]*ys[i]*ys[i];
	}
    double sigma_x2y2 = arr_sum(tmp,n);
    
  
    //sigma_z1x1
    for(int i=0;i<n;i++)
    {
    	tmp[i] = zs[i]*xs[i];
	}
    double sigma_z1x1= arr_sum(tmp,n);    
    
 
    //sigma_z1y1
    for(int i=0;i<n;i++)
    {
    	tmp[i] = zs[i]*ys[i];
	}
    double sigma_z1y1= arr_sum(tmp,n); 
    
    
    //sigma_z1x2
    for(int i=0;i<n;i++)
    {
    	tmp[i] = zs[i]*xs[i]*xs[i];
	}
    double sigma_z1x2= arr_sum(tmp,n);     
    

    //sigma_z1y2
    for(int i=0;i<n;i++)
    {
    	tmp[i] = zs[i]*ys[i]*ys[i];
	}
    double sigma_z1y2= arr_sum(tmp,n);         
    
    
    //sigma_z1x1y1
    for(int i=0;i<n;i++)
    {
    	tmp[i] = zs[i]*xs[i]*ys[i];
	}
    double sigma_z1x1y1= arr_sum(tmp,n);      
    
    
    double a[6 * 6]={
			sigma_x4,sigma_x2y2,sigma_x3y1,sigma_x3,sigma_x2y1,sigma_x2,
			sigma_x2y2,sigma_y4,sigma_x1y3,sigma_x1y2,sigma_y4,sigma_y2,
			sigma_x3y1,sigma_x1y3,sigma_x2y2,sigma_x2y1,sigma_x1y2,sigma_x1y1,
			sigma_x3,sigma_x1y2,sigma_x2y1,sigma_x2,sigma_x1y1,sigma_x1,
			sigma_x2y1,sigma_y4,sigma_x1y2,sigma_x1y1,sigma_y2,sigma_y1,
			sigma_x2,sigma_y2,sigma_x1y1,sigma_x1,sigma_y1,5,
			};
	    
    double b[6]  = {sigma_z1x2,sigma_z1y2,sigma_z1x1y1,sigma_z1x1,sigma_z1y1,sigma_z1};

    double t[6]  = {};
    
    gauss(a, b, t, 6);
    
	/*
    printf("增广矩阵\n");
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            printf("%f\t", a[i * 6 + j]);
        }
        printf("%f\n", b[i]);
    }
    
    printf("\n");
    
    printf("解向量\n");
    for (int i = 0; i < 6; i++)
        printf("%f\t", t[i]);
    printf("\n\n");
    */
    
    //生成粗糙曲面
	double p[256] = {0};
	double q[256] = {0};  
    for(int i=0;i<256;i++)
    {
    	p[i] = q[i] = -20+0.15625*i;
	}
    for(int i=0;i<256;i++)
    {
    	for(int j=0;j<256;j++)
    	{	
    		double x = p[i];
    		double y = q[j];
    		double z = x*x*t[0]+y*y*(0.99*t[0]+t[1])+x*y*t[2]+x*t[3]+y*t[4]+t[5];
    		surface[i*256+j] = 0.1*z*sin(0.01*x*y)+0.6*exp(sin(0.6*cos(x+y)))-0.3*exp(sin(0.3*exp(sin(x))-3*sin(y)))+z;
		}
	}
    
    //加入凹凸算子
	for(int num=0;num<100;num++)
	{
    	int size  = random(10,30);
    	double scale = random(2,8)/10;
    	double patch[size*size] = {0};
    	for(int i=0;i<size*size;i++)
		{
			patch[i] = 1;
		}
		//生成凹凸块 
		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++)
			{
				int idx = i*size + j;
				for(int k=0;k<(size+1)/2;k++)
				{	
					if(i==k||i==(size-k)||j==k||j==(size-k))
					{
						patch[idx] = 1+0.1*k;
						break;
					}
				}
			}
		
		}
		
    	/*
		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++)
			{
				printf(" %f ",patch[i*size+j]);
			}
			printf("\n\n");
		}
	     */ 
	     
    	int start_x = random(0,255-size);
    	int start_y = random(0,255-size);
    	/*
    	printf("start_x:%d\n",start_x);
    	printf("start_y:%d\n",start_y);
    	*/
    	//加入凹凸块 
    	for(int i=0;i<size;i++)
    	{
    		for(int j=0;j<size;j++)
    		{
    			int patch_idx = i*size + j;
    			int surface_idx = (start_y+i)*256+start_x+j;
    			surface[surface_idx] += scale*patch[patch_idx];
			}
		}
    	
    	//printf("\n\n");
	}
}

double random(int min,int max)
{	
	//srand((unsigned)time(NULL));
	return rand()%(max-min+1)+min;
}

double arr_sum(double xs[],int n)
{	
	double res = 0;
	for(int i=0;i<n;i++)
	{
		res += xs[i];
	}
	return res;
}

void gauss(double a[], double b[], double x[], int n)//列主元高斯消去法，a[]系数矩阵，b[]结果向量，x[]解向量，n阶数
{
    int i, j, k, exchangeline, exchangeflag = 0;
    double temp, max;

    for (k = 0; k < n - 1; k++) { //k迭代次数
        max = a[k * n + k];
        for (i = k + 1; i < n; i++) { //寻找主元，i行号
            if (fabs(max) < fabs(a[n * i + k])) {
                max          = a[n * i + k];
                exchangeflag = 1;     //交换标志
                exchangeline = n * i; //记录需要交换的行号
            }
        }
        if (exchangeflag) { //换行，j列号
            for (j = 0; j < n; j++) {
                temp                = a[exchangeline + j]; //对系数矩阵操作
                a[exchangeline + j] = a[k * n + j];https://github.com/Reyyyyyy/Simulation-surface-generating-program-based-on-C-and-Python/blob/main/algorithm.cpp
                a[k * n + j]        = temp;
            }
            temp                = b[exchangeline / n]; //对结果向量操作
            b[exchangeline / n] = b[k];
            b[k]                = temp;
            exchangeflag        = 0; //清除交换标志
        }
        for (i = k + 1; i < n; i++) { //消元
            temp = a[i * n + k] / a[k * n + k];
            b[i] = b[i] - b[k] * temp; //对结果向量操作
            for (j = k; j < n; j++)
                a[i * n + j] = a[i * n + j] - a[k * n + j] * temp; //对系数矩阵操作
        }
    }
    for (i = n - 1; i > -1; i--) { //回代
        temp = b[i];
        for (j = n - 1; j > i; j--)
            temp = temp - a[n * i + j] * x[j];
        x[i] = temp / a[i * n + i];
    }
}



