#pragma once

//种群大小
const int popsize = 100;
//编码长度  一次发车的数量
int chromlength = 10;
//候选数量
int candidate = 2;
//交叉概率
const float pc = 0.3;
//变异概率
const float pm = 0.001;
//最大迭代次数
int iteration = 40;
//提前结束循环
int early_stop = 10;


//初始化种群
int* initpop(int popsize,int chromlength);
//Genetic algorithms 遗传算法主循环函数
void GAmain(double cr0);                  
//Mutation operation 变异操作
int * mutation(int * pop, float pm);		
//Crossover Operation 交叉操作
int * crossover(int * pop, float pc);		
//选择复制操作
float* selection(float* pop, float* fitvalue); 


double myu(double a, double b) // Uniform Distribution
{
	double y;
	if (a>b) {
		printf("\nThe first parameter should be less than the second!");
		exit(1);
	}
	y = (double)rand() / (RAND_MAX);
	return (a + (b - a)*y);
}

//产生a-b之间的随机数
int rdint(int a, int b)
{
	double x;
	int temp;
	x = myu(a, b + 0.99);
	temp = int(floor(x));
	return(temp);
}