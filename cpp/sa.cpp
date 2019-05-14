#include<iostream>
#include<algorithm>
#include<vector>
#include<numeric>
#include<math.h>
#include<time.h>
#include"sa.h"

float* binary2decimal(int* pop)
{
	float * newpop = new float[popsize * chromlength];
	for (int i = 0; i < chromlength; i++) {
		for (int j = 0; j < popsize; j++) {
			newpop[j * chromlength + i] = pow((float)2.0, chromlength - i - 1) * (float)pop[j * chromlength + i];
		}
	}
	float * res = new float[popsize];
	memset(res, 0, sizeof(float)*popsize);
	for (int i = 0; i < chromlength; i++) {
		for (int j = 0; j < popsize; j++) {
			res[j] = res[j] + newpop[j * chromlength + i];
		}
	}
	for (int j = 0; j < popsize; j++) {
		res[j] = res[j] * 10.0 / 1023.0;
	}
	return res;
}

float binary2decimalsingle(int* pop)
{
	float * newpop = new float[chromlength];
	for (int i = 0; i < chromlength; i++) {
		newpop[i] = pow((float)2.0, chromlength - i) * (float)pop[i];
	}
	float res = 0.0;
	for (int i = 0; i < chromlength; i++) {
		res = res + newpop[i];
	}
	res = res * 10.0 / 1023.0;

	return res;
}

float * cal_objvalue(int * pop)
{
	float * objvalue = new float[popsize];
	float * x = binary2decimal(pop);
	/*for (int j = 0; j < popsize; j++) {
		cout << x[j] << " ";
	}
	cout << endl;*/
	for (int j = 0; j < popsize; j++) {
		objvalue[j] = 10.0 * sin(5.0 * x[j]) + 7.0 * abs(x[j] - 5) + 10.0;
		//cout << x[j] << " " << objvalue[j] << endl;

	}
	//double temp = 10.0 * sin(5.0 * 4) + 7.0 * abs(4 - 5) + 10.0;
	//cout << temp << endl;
	return objvalue;
}

///初始化种群
int* initpop(int popsize, int chromlength)
{
	//float smoth = 10;
	int* pop = new int[popsize  * chromlength];
	for (int i = 0; i < popsize; i++) {
		for (int j = 0; j < chromlength; j++) {
			pop[i * chromlength + j] = rdint(0,1);
		}
	}
	return pop;
}
//选择复制操作
//fitvalue  大小：popsize * 1  种群大小
int* selection(int* pop, float* fitvalue)
{
	int * newpop = new int[popsize * chromlength];
	float totalfit = 0;
	for (int i = 0; i < popsize; i++) {
		totalfit += fitvalue[i];
	}
	float * new_fitvalue = new float[popsize];
	float * ms = new float[popsize];
	for (int i = 0; i < popsize; i++) {
		if (i == 0) {
			new_fitvalue[i] = fitvalue[i] / totalfit;
		}
		else {
			new_fitvalue[i] = new_fitvalue[i - 1] + fitvalue[i] / totalfit;
		}
		ms[i] = float(rand() % 10000) / 10000.0;
	}
	std::sort(ms, ms + popsize);
	int fitin = 0, newin = 0;
	while (newin < popsize) {
		if (ms[newin] < new_fitvalue[fitin]) {
			for (int j = 0; j < chromlength; j++) {
				newpop[newin * chromlength + j] = pop[fitin * chromlength + j];
			}
			newin = newin + 1;
		}
		else {
			fitin = fitin + 1;
		}
	}
	return newpop;
}

//Crossover Operation 交叉操作
int * crossover(int * pop, float pc)
{
	int * newpop = new int[popsize * chromlength];
	for (int i = 0; i < popsize - 1; i += 2) {
		float rand01 = (rand() % 10000) / 10000.0;
		if (rand01 < pc) {
			int cpoint = round((rand() % 10000) / 10000.0 * chromlength);
			//核对
			for (int j = 0; j < cpoint; j++) {
				newpop[i * chromlength + j] = pop[(i + 1) * chromlength + j];
				newpop[(i + 1) * chromlength + j] = pop[i * chromlength + j];
			}
			for (int j = cpoint; j < chromlength; j++) {
				newpop[i * chromlength + j] = pop[i * chromlength + j];
				newpop[(i + 1) * chromlength + j] = pop[(i + 1) * chromlength + j];
			}
		}
		else {
			for (int j = 0; j < chromlength; j++) {
				newpop[i * chromlength + j] = pop[i * chromlength + j];
				newpop[(i + 1) * chromlength + j] = pop[(i + 1) * chromlength + j];
			}
		}
	}
	return newpop;
}

//Mutation operation 变异操作
int * mutation(int * pop, float pm)
{
	int * newpop = new int[popsize * chromlength];
	for (int i = 0; i < popsize; i++) {
		float rand01 = (rand() % 10000) / 10000.0;
		if (rand01 < pm) {
			for (int j = 0; j < chromlength; j++) {
				newpop[i * chromlength + j] = pop[i * chromlength + j];
			}
			int mpoint = rand() % chromlength;		//变异的位置
			int mutation_num = rand() % candidate;	//变异后的数值
			newpop[i * chromlength + mpoint] = mutation_num;
		}
		else {
			for (int j = 0; j < chromlength; j++) {
				newpop[i * chromlength + j] = pop[i * chromlength + j];
			}
		}
	}

	return newpop;
}

//找最优解
int * get_best(int * pop, float * fitvalue, float &res)
{
	float bestfit = fitvalue[0];
	int * bestindividual = new int[chromlength];
	for (int i = 0; i < chromlength; i++) {
		bestindividual[i] = pop[i];
	}
	for (int i = 1; i < popsize; i++) {
		if (fitvalue[i] > bestfit) {
			for (int j = 0; j < chromlength; j++) {
				bestindividual[j] = pop[i * chromlength + j];
			}
			bestfit = fitvalue[i];
		}
	}
	res = bestfit;
	return bestindividual;
}

void GAmain()       //Genetic algorithms 遗传算法
{
	srand((unsigned)time(NULL));
	int* pop = initpop(popsize, chromlength);
	int i = iteration;

	float bestfit = 0;
	int * best_answer = NULL;
	int patience = early_stop;
	while (i--)
	{
		//计算适应度值（函数值）
		float * fitvalue = cal_objvalue(pop);

		//选择进行复制 轮盘赌方式
		int * newpop = selection(pop, fitvalue);

		//交叉操作
		newpop = crossover(newpop, pc);

		//变异操作
		newpop = mutation(newpop, pm);

		//更新种群
		pop = newpop;

		//寻找最优解，并记录下来
		float curr_bestfit = 0;
		int * res = get_best(pop, fitvalue, curr_bestfit);

		if (curr_bestfit > bestfit) {
			bestfit = curr_bestfit;
			best_answer = res;
			patience = early_stop;
		} else {
			patience -= 1;
		}

		std::cout << "current best value is:" << bestfit << std::endl;
		
		//如果连续 early_stop 次找不到更优解，停止迭代，节约时间
		if (patience == 0) {
			std::cout << "early_stop,best value is:" << bestfit << std::endl;
			break;
		}
	}
}

int main()
{
	GAmain();
	return 0;
}