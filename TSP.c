#include <stdio.h>
#define _CRT_SECURE_NO_WARNINGS
#include "mpi.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define NPROC 17     //进程数量（主从各一个宏定义）@@@@@@@@@@@@@
#define N_COLONY 100 // N_COLONY>=xColony300
#define CITY 783
//#define TWOGRP 20               //组间交换和漫游比例是19：1这个问题值得考虑修改？！！实际是20：1！！
#define GRPSCALE 4                //组的规模@@@@@@@@@@@@
#define grpnum (NPROC / GRPSCALE) //组的个数

#define GLOBAL_EXPERIMENT_TIME 30 //重复实验次数

double sumTemp, sumbest;
int xColony = N_COLONY;      //##//  个体数
int xCity = CITY;            //（问题相关）
double edgeSpeed = 5000;     //##//  临界速度
double const PROBAB2 = 0.04; //##//  映射概率 //0.04(80) 0.03(50) 0.015(80)  0.05(50);
#define DIVISION 4
const double PROBAB1 = PROBAB2 / DIVISION; //##//  变异概率

const double PROBAB3 = 0.001;   //##//  映射变异临界值
const double PROBAB4 = 0.00256; //##//  聚类评价函数变化临界值

long Ni, NOCHANGE = 2000000; //##//  最大停止改变代数
long loopcounter = 0;
long MAXGEN = 10000000; //5
long INTERVAL = 5000;
int colony[N_COLONY][CITY]; //所有个体的染色体
double cityXY[CITY][2];
double city_dis[CITY][CITY];
double dis_p[N_COLONY]; //所有个体的评估值
double tempdis;         //传回的个体的评估值
int tempcol[CITY];      //传回的个体的染色体
double speed;
int temp[CITY], ibest, iworst, ipass;
clock_t timeStart, timeNow, timeTemp;
char filepath[100], filepath2[100], filepath3[100];

void initm();
void init();

int position(int *tmp, int C);
void invert(int pos_start, int pos_end, int *ptr = temp);
void reverse(int *ptr, int pos);
void topo();
void printBest(int g_looptimer, int GenNum, int Ni);
void tempTest(int i, int *ptr = temp);
void mapped();
//void LastCP();
double path(int tmp[], int k1, int k2);
double SPAD_compute();
int curBest[CITY];
int **mat = NULL;
int **matlocalbest = NULL;

static union
{
    double *best_ptr = (double *)0x7f7fffffffffffff;
    double best; //初始化为最大值
};

int main(int argc, char *argv[])
{
    int r;
    int mytid, numprocs;
    char mname[30];
    MPI_Init(&argc, &argv); /**/                /*程序初始化*/
    MPI_Comm_rank(MPI_COMM_WORLD, &mytid); /**/ /*得到当前进程号*/
    MPI_Get_processor_name(mname, &r);
    mat = (int **)malloc(sizeof(int *) * CITY);
    matlocalbest = (int **)malloc(sizeof(int *) * CITY);
    for (int i = 0; i < CITY; ++i)
    {
        mat[i] = (int *)malloc(sizeof(int) * CITY);
        matlocalbest[i] = (int *)malloc(sizeof(int) * CITY);
    }
    MPI_Status status;
    int sign = 1, i, i1, j1;
    double result;
    double data[NPROC], tempdata, *ptrdata[NPROC], dt;       //指针数组用来存放精英的顺序
    int chrom[NPROC][CITY], tempchr[CITY], *ptrchrom[NPROC]; //（问题相关）第一维是进程数，第二维比城市数大。
    int j, k, m, x;                                          //用于交换，m用来做记录组别的变量
int times = 0;                                           //times是提交的次数

    FILE *fpme;
    double t, probab1 = PROBAB1, probab2 = PROBAB2;
    register int C1, js, ks, pos_C, pos_C1;
    int k1, k2, l1, l2, pos_flag;
    register double disChange;
    static int is = 0;
    FILE *fppp;
    double SPAD0 = 1.0, SPAD = 1.0; //,aveSPAD;
    srand((unsigned)time(NULL) + mytid);
    strcpy(filepath, "./ALL_tsp/rat783.tsp");
    //strcpy_s(filepath, argv[1]);
    if (mytid == 0)
    {
        //strcpy(filepath2, argv[2]);
        strcpy(filepath2, "./result2.txt");
        strcpy(filepath3, "./result1.txt");
        //初始化指针数组
        for (i = 0; i < NPROC; i++)
        {
            ptrdata[i] = &data[i];
            ptrchrom[i] = chrom[i];
        }
        //初始化指针数组
        initm();
        for (i = 1; i < NPROC; i++)
            MPI_Send(cityXY, CITY * 2, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(cityXY, CITY * 2, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
    }
    for (int g_looptimer = 0; g_looptimer < GLOBAL_EXPERIMENT_TIME; ++g_looptimer)
    {
        init();
        best_ptr = (double *)0x7f7fffffffffffff; // 等价于将best设置为最大
        times = 0;
        //以上初始化
        for (loopcounter = 0; loopcounter <= MAXGEN; loopcounter++)
        {
            if (loopcounter % INTERVAL == 0 && loopcounter != 0)
            {
                times++;
                //第二次交流////////////////////////
                //收
                if (mytid == 0)
                {
                    for (i = 1; i < NPROC; i++) //接收最优染色体为当前解
                    {
                        MPI_Recv(&data[i - 1], 1, MPI_DOUBLE, i, 97, MPI_COMM_WORLD, &status);
                        MPI_Recv(chrom[i - 1], CITY, MPI_INT, i, 97, MPI_COMM_WORLD, &status);
                        //printf("I got %f from %d\n",data[i],i+1);
                    }
                    for (i = 0; i < NPROC - 1; i++)
                    {
                        if (data[i] < best)
                        {
                            best = data[i];
                            for (k = 0; k < xCity; ++k)
                            {
                                curBest[k] = chrom[i][k];
                            }
                        }
                        printf("from colony%d %lf\n", i + 1, data[i]);
                    }
                    printf("Time(s)=%d,epoch%d,best-so-far is %lf\n", times, loopcounter, best);
                    printf("Cluster Evaluate %lf->%lf\n", SPAD0, SPAD);
                    for (i = 1; i < NPROC; i++)
                    {
                        MPI_Recv(&data[i - 1], 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status);
                        MPI_Recv(chrom[i - 1], CITY, MPI_INT, i, 99, MPI_COMM_WORLD, &status);
                        //printf("I got %f from %d\n",data[i],i+1);
                    }
                    int group_count = grpnum;
                    if (fabs(SPAD - SPAD0) >= PROBAB4)
                    {
                        for (i = 0; i < NPROC - 1; i++) //1-NPROC
                        {
                            int m = i / group_count;
                            m = m * GRPSCALE + ((i + 1) % GRPSCALE) + 1;
                            MPI_Send(&data[i], 1, MPI_DOUBLE, m, 99, MPI_COMM_WORLD);
                            MPI_Send(chrom[i], CITY, MPI_INT, m, 99, MPI_COMM_WORLD);
                            //printf("I send %f to %d\n",data[i],i+1);
                        }
                    }
                    else
                    {
                        for (i = 0; i < NPROC - 1; i++) //1-NPROC
                        {
                            int m = i / group_count;
                            m = ((m + 1) * GRPSCALE + ((i + 1) % GRPSCALE)) % (NPROC - 1) + 1; //群落j的第i个种群向群落j+1的i+1种群迁移
                            MPI_Send(&data[i], 1, MPI_DOUBLE, m, 99, MPI_COMM_WORLD);
                            MPI_Send(chrom[i], CITY, MPI_INT, m, 99, MPI_COMM_WORLD);
                            //printf("I send %f to %d\n",data[i],i+1);
                        }
                    }
                    SPAD0 = SPAD;
                    double SPAD_SET[NPROC];
                    for (i = 1; i < NPROC; i++)
                    {
                        MPI_Recv(SPAD_SET + i - 1, 1, MPI_DOUBLE, i, 98, MPI_COMM_WORLD, &status);
                    }
                    double aver = 0, group_aver[grpnum] = {0.0}, a1 = 0, a2 = 0, tmp;
                    for (j = 0; j < group_count; ++j)
                    {
                        for (i = 0; i < GRPSCALE; ++i)
                        {
                            group_aver[j] += SPAD_SET[j * GRPSCALE + i];
                        }
                        aver += group_aver[j];
                        group_aver[j] /= GRPSCALE;
                    }
                    aver /= (NPROC - 1);
                    for (j = 0; j < group_count; ++j)
                    {
                        tmp = (group_aver[j] - aver);
                        a1 += tmp * tmp;
                    }
                    for (i = 0; i < NPROC - 1; ++i)
                    {
                        tmp = (SPAD_SET[i] - aver);
                        a2 += tmp * tmp;
                    }
                    SPAD = fabs(a2) < 1e-9 ? 1.0 : a1 * (NPROC - 1 - group_count) / (a2 * (NPROC - 2));
                }
                else
                {
                    //发最优染色体给主进程0
                    SPAD0 = SPAD;
                    SPAD = SPAD_compute();
                    ipass = ibest;
                    MPI_Send(&dis_p[ipass], 1, MPI_DOUBLE, 0, 97, MPI_COMM_WORLD);
                    MPI_Send(colony[ipass], CITY, MPI_INT, 0, 97, MPI_COMM_WORLD);
                    //发一条随机染色体给主进程0
                    ipass = rand() % xColony;
                    MPI_Send(&dis_p[ipass], 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
                    MPI_Send(colony[ipass], CITY, MPI_INT, 0, 99, MPI_COMM_WORLD);
                    //收@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    //发
                    //从主进程要回群落内另一个种群的染色体
                    MPI_Recv(&tempdis, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
                    MPI_Recv(tempcol, CITY, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
                    MPI_Send(&SPAD, 1, MPI_DOUBLE, 0, 98, MPI_COMM_WORLD);
                    double t = DIVISION * SPAD;
                    probab1 = PROBAB1 * (1 - t + t * t / 2 - t * t * t / 6); //泰勒展开取近似
                    // *e^(-SPAD*division) 内逆转概率随种群内相似度升高而升高
                    t -= DIVISION;
                    probab2 = PROBAB2 * (1 + t + t * t / 2 + t * t * t / 6);
                    // *e^(division*(SPAD-1)) 部分交换概率随种群内差异度升高而升高
                    if (rand() % 1000 / 1000 < (1 - SPAD) || fabs(SPAD - SPAD0) < 1e-9) //迁移
                    {
                        ipass = rand() % xColony;
                        dis_p[ipass] = tempdis;
                        for (j1 = 0; j1 < xCity; j1++)
                            colony[ipass][j1] = tempcol[j1];
                    }
                    else
                    {
                        int ranStart = rand() % xCity, pos;
                        for (i = 0; i < xColony; ++i)
                        {
                            pos = position(colony[i], ranStart);
                            reverse(colony[i], pos); // 把pos处的城市通过翻转拉到染色体头部
                        }
                    }
                }
            }
            if (mytid != 0) //子进程程序
            {
                ///第二次交流
                ////////////////////////////
                for (js = 0; js < xCity; js++)
                    temp[js] = colony[is][js];
                disChange = 0;
                pos_flag = 0;
                pos_C = rand() % xCity;
                ////////////////////////////
                for (;;)
                {
                    //两种算子至少进行一个， 应该除以两种变异至少发生一个的概率
                    if ((rand() * 1.0 / RAND_MAX) * (1 - (1 - probab1) * (1 - probab2)) < probab1)
                    {
                        do
                            pos_C1 = rand() % xCity;
                        while (pos_C1 == pos_C);
                        C1 = colony[is][pos_C1];
                    } //if
                    else
                    {
                        do
                            js = rand() % xColony;
                        while (js == is);
                        ks = position(colony[js], temp[pos_C]);
                        C1 = colony[js][(ks + 1) % xCity];
                        pos_C1 = position(temp, C1);
                    }
                    if (speed > edgeSpeed && pos_C1 < pos_C + 2)
                        break; ///////////////////////
                    if ((pos_C + 1) % xCity == pos_C1 || (pos_C - 1 + xCity) % xCity == pos_C1)
                        break;
                    k1 = temp[pos_C];
                    k2 = temp[(pos_C + 1) % xCity];
                    l1 = temp[pos_C1];
                    l2 = temp[(pos_C1 + 1) % xCity];
                    disChange += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
                    invert(pos_C, pos_C1);
                    pos_flag++;
                    if (pos_flag > xCity - 1)
                        break; ////////////
                    pos_C++;
                    if (pos_C >= xCity)
                        pos_C = 0; /**********************/
                    if (speed < edgeSpeed && disChange < 0)
                    {
                        dis_p[is] += disChange;
                        disChange = 0;
                        tempTest(is);
                    } //每有改变就测试是否是目前最优解
                }
                if (speed >= edgeSpeed && disChange < 0)
                {
                    dis_p[is] += disChange;
                    disChange = 0;
                    tempTest(is);
                } /////speed>=1500 &&
                is++;
                if (is >= xColony)
                {
                    Ni++;
                    is = 0;
                    if ((rand() * 1.0 / RAND_MAX < probab2)) // speed<edgeSpeed &&
                    {
                        if ((rand() & 32767) / 32767.0 < PROBAB3 + fabs(SPAD - SPAD0))
                            mapped(); //映射算子
                        topo();       //拓扑算子
                    }
                }
            } //while(loopcounter++<=MAXGEN);
        }
        if (mytid == 0)
        {
            fpme = fopen(filepath2, "a");
            printf("This is a result of %d:%lf\n\n\n", CITY, best);
            fprintf(fpme, "%lf\n", CITY, best);
            fclose(fpme);
            printBest(g_looptimer, loopcounter, Ni);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    for (i = 0; i < CITY; i++)
    {
        free(mat[i]);
        free(matlocalbest[i]);
    }
    free(mat);
    free(matlocalbest);
    MPI_Finalize();
    return 0;
}
double SPAD_compute() //对边记数，计算不同染色体群落中边集与该群落中最优的染色体的边集之间的差异
{
    int ibest2 = 0;
    int xx, yy, zz;
    //int mat[CITY][CITY], matlocalbest[CITY][CITY];
    int best_i, count_of_1;
    double D[N_COLONY], STD, lSTD;
    for (xx = 0; xx < N_COLONY; xx++)
    {
        if (dis_p[ibest2] > dis_p[xx])
            ibest2 = xx;
    }
    ibest = ibest2;
    for (xx = 0; xx < xCity; xx++)
        for (yy = 0; yy < xCity; yy++)
            matlocalbest[xx][yy] = 0;
    for (xx = 0; xx < xCity; xx++)
    {
        matlocalbest[colony[ibest2][xx]][colony[ibest2][(xx + 1) % xCity]] = 1;
    }
    lSTD = 0;
    for (zz = 0; zz < N_COLONY; zz++)
    {
        for (xx = 0; xx < CITY; xx++) //ci chu gui 0
        {
            for (yy = 0; yy < CITY; yy++)
                mat[xx][yy] = 0;
        }
        for (xx = 0; xx < xCity; xx++)
        {
            mat[colony[zz][xx]][colony[zz][(xx + 1) % xCity]] = 1; //全部解矩阵建立
        }
        count_of_1 = 0;
        for (xx = 0; xx < xCity; xx++)
        {
            for (yy = 0; yy < xCity; yy++)
            {
                if (matlocalbest[xx][yy] == 1 && mat[xx][yy] == 1)
                {
                    break;
                }
            }
            if (yy < xCity)
                count_of_1++;
        }
        //求出了一个个体与本地最好个体的差异
        lSTD += xCity - count_of_1; //累加差异
    }                               //for zz
    return lSTD / (N_COLONY * xCity);
}
void reverse(int *ptr, int pos) //翻转染色体以使pos对齐开头
{
    int j, k, t;
    if (pos <= 0)
        return;
    j = 0;
    k = pos;
    for (; j < k; j++, k--)
    {
        t = ptr[j];
        ptr[j] = ptr[k];
        ptr[k] = t;
    }
    if (pos < xCity - 1)
    {
        j = pos + 1;
        k = xCity - 1;
        for (; j < k; j++, k--)
        {
            t = ptr[j];
            ptr[j] = ptr[k];
            ptr[k] = t;
        }
    }
}
void invert(int pos_start, int pos_end, int *ptr)
{
    int j, k, t;
    if (pos_start < pos_end)
    {
        j = pos_start + 1;
        k = pos_end;
        for (; j <= k; j++, k--)
        {
            t = ptr[j];
            ptr[j] = ptr[k];
            ptr[k] = t;
        }
    }
    else
    {
        if (xCity - 1 - pos_start <= pos_end + 1)
        {
            j = pos_end;
            k = pos_start + 1;
            for (; k < xCity; j--, k++)
            {
                t = ptr[j];
                ptr[j] = ptr[k];
                ptr[k] = t;
            }
            k = 0;
            for (; k <= j; k++, j--)
            {
                t = ptr[j];
                ptr[j] = ptr[k];
                ptr[k] = t;
            }
        }
        else
        {
            j = pos_end;
            k = pos_start + 1;
            for (; j >= 0; j--, k++)
            {
                t = ptr[j];
                ptr[j] = ptr[k];
                ptr[k] = t;
            }
            j = xCity - 1;
            for (; k <= j; k++, j--)
            {
                t = ptr[j];
                ptr[j] = ptr[k];
                ptr[k] = t;
            }
        }
    }
}
int position(int *tmp, int C)
{
    int j;
    for (j = 0; j < xCity; j++)
        if (*(tmp + j) == C)
            break;
    return (j);
}
void tempTest(int i, int *ptr)
{
    int j;
    double dt;
    if (ptr == temp)
        for (j = 0; j < xCity; j++)
            colony[i][j] = ptr[j];
    if ((int)sumbest > (int)dis_p[i])
    {
        sumbest = dis_p[i];
        ibest = i;
        Ni = 0;
        timeNow = clock();
        dt = (double)(timeNow - timeTemp) / CLOCKS_PER_SEC;

        if (dt > 0.1)
        {
            speed = (sumTemp - sumbest) / dt;
            sumTemp = sumbest;
            timeTemp = timeNow;
        }
        //printf("\n%f   %4.2f  .1f",sumbest,(double)(timeNow-timeStart)/CLOCKS_PER_SEC,speed);
    }
}
void printBest(int timer, int GenNum, int Ni)
{
    FILE *fpme = fopen(filepath3, "a");
    fprintf(fpme, "Experiment %d : %lf\n", timer + 1, best);
    fprintf(fpme, "\n   CITY      %d\t\tN_COLONY  %d", CITY, N_COLONY);
    fprintf(fpme, "\n   maxGen    %d\t\ttime      %4.2f  seconds", MAXGEN, (double)(timeNow - timeStart) / CLOCKS_PER_SEC);
    fprintf(fpme, "\n   GenNum    %d\t\tNo change %Ld\n\n", GenNum, Ni);
    for (int i = 0; i < xCity; i++)
    {
        if (i % 10 == 0 && i != 0)
        {
            fprintf(fpme, "\n");
        }
        fprintf(fpme, "]->", curBest[i] + 1);
    }
    fprintf(fpme, "\n\n");
    fclose(fpme);
}
void mapped()
{
    int start, end, i, j, k, kt, t, disPlace, kDC, kC;
    i = rand() % xColony;
    j = rand() % xColony;
    if (i == j)
        return;
    if (dis_p[i] < dis_p[j])
    {
        t = i;
        i = j;
        j = t;
    }
    for (k = 0; k < xCity; k++)
        temp[k] = colony[i][k];
    /////////////////
    start = rand() % xCity;
    end = (start + rand() % 180 + 20) % xCity; //rand()%xCity;
    kt = position(temp, colony[j][start]);     //部分映射一二位同
    disPlace = kt - start;
    if (temp[(kt + 1) % xCity] == colony[j][(start + 1) % xCity])
    {
        if (start >= end)
            end += xCity;
        for (k = start; k <= xCity; k++)
        {
            kDC = (k + disPlace) % xCity;
            kC = k % xCity;
            if (temp[kDC] == colony[j][kC])
                continue;
            t = position(temp, colony[j][kC]);
            temp[t] = temp[kDC];
            temp[kDC] = colony[j][kC];
        }
    }
    else
    {
        if (temp[(kt - 1 + xCity) % xCity] == colony[j][(start + 1) % xCity])
        {
            if (start >= end)
                end += xCity;
            for (k = kt = start; k <= end; k++, kt--)
            {
                kDC = (kt + xCity + disPlace) % xCity;
                kC = k % xCity;
                if (temp[kDC] == colony[j][kC])
                    continue;
                t = position(temp, colony[j][kC]);
                temp[t] = temp[kDC];
                temp[kDC] = colony[j][kC];
            }
        }
        else
            return;
    }
    double temp_dis = 0;
    for (j = 0; j < xCity - 1; j++)
        temp_dis += city_dis[temp[j]][temp[j + 1]];
    temp_dis += city_dis[temp[0]][temp[xCity - 1]];
    dis_p[i] = temp_dis;
    /*********/
    tempTest(i);
}
void topo()
{
    int i, j, k, t;
    i = rand() % xColony;
    do
        j = rand() % xColony;
    while (i == j);
    if (dis_p[i] < dis_p[j])
    {
        t = i;
        i = j;
        j = t;
    }
    int index[CITY], shift[CITY], indexCount[CITY] = {0}, Circle[CITY] = {-1}, circleCount = 0;
    for (k = 0; k < xCity; ++k)
        index[colony[i][k]] = k;
    for (k = 0; k < xCity; ++k)
    {
        t = k;
        if (Circle[k] == -1) //通过倒索引建立拓扑环
        {
            do
            {
                Circle[t] = circleCount;
                t = index[colony[j][t]];
            } while (t = k);
            ++circleCount;
        }
    }
    int key1, key2;
    int oppo_key1, oppo_key2;
    int k1, k2, l1, l2;
    double disChange1, disChange2;
    if (circleCount > 1) // 环数大于1时才能拓扑合并
    {
        int cir2 = 1 + (rand() % (circleCount - 1)), cir1 = cir2 >> 1; // 选择两个环合并
        key1 = -1;
        key2 = -1;
        for (k = 0; k < xCity; ++k)
        {
            if (Circle[k] == cir2)
            {
                key2 = k;
                break;
            }
        }
        for (k = xCity - 1; k >= 0; --k)
        {
            if (Circle[k] == cir1)
            {
                key1 = k;
                break;
            }
        }
        oppo_key1 = index[colony[j][key1]];
        oppo_key2 = index[colony[j][key2]];
        k1 = colony[j][(key1 - 1 + xCity) % xCity];
        k2 = colony[j][key1];
        l1 = colony[j][key2];
        l2 = colony[j][(key2 + 1) % xCity];
        disChange1 = city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
        k1 = colony[i][(oppo_key1 - 1 + xCity) % xCity];
        k2 = colony[i][oppo_key1];
        l1 = colony[i][oppo_key2];
        l2 = colony[i][(oppo_key2 + 1) % xCity];
        disChange1 = city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
        if (disChange1 < 0)
        {
            dis_p[j] += disChange1;
            invert(key1 - 1, key2, colony[j]);
        }
        if (disChange2 < 0)
        {
            dis_p[i] += disChange2;
            invert(oppo_key1 - 1, oppo_key2, colony[i]);
            index[colony[j][key1]] = oppo_key2;
            index[colony[j][key2]] = oppo_key1;
        }
    }
    ///注意！ 拓扑合并后Circle的值不保证正确！
    ///变异翻转
    for (k = 0; k < xCity; ++k)
    {
        t = shift[k] = (index[colony[j][k]] - k + xCity) % xCity; //计算每一位的偏移
        ++indexCount[t];
    }
    int most = 0, _max = 0;
    bool flag = false;
    for (k = 0; k < xCity; ++k)
    {
        if (indexCount[k] > _max)
        {
            most = k;
            _max = indexCount[k]; //获取出现次数最多的偏移量
        }
    }
    int l = 0, r = -1;
    for (k = 0; k < xCity; ++k)
    {
        if (shift[k] == most)
        {
            flag = true;
            if (k == r + 1)
            {
                ++r;
                if (r - l + 1 >= (xCity >> 4)) // 不相同区间超过了染色体的1/16 直接选择倒转
                    break;
            }
            else
            {
                if (r - l > 0 && rand() % xCity < (r - l + 1))
                    break;
                l = r = k;
            }
        }
        else if (flag)
            break;
    }
    if (r - l <= 0)
        return;
    key1 = l;
    key2 = r;
    oppo_key1 = index[colony[j][key1]];
    oppo_key2 = index[colony[j][key2]];
    k1 = colony[j][(key1 - 1 + xCity) % xCity];
    k2 = colony[j][key1];
    l1 = colony[j][key2];
    l2 = colony[j][(key2 + 1) % xCity];
    disChange1 = city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
    k1 = colony[i][(oppo_key1 - 1 + xCity) % xCity];
    k2 = colony[i][oppo_key1];
    l1 = colony[i][oppo_key2];
    l2 = colony[i][(oppo_key2 + 1) % xCity];
    disChange2 = city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
    if (disChange1 < 0)
    {
        dis_p[j] += disChange1;
        invert(key1 - 1, key2, colony[j]);
    }
    if (disChange2 < 0)
    {
        dis_p[i] += disChange2;
        invert(oppo_key1 - 1, oppo_key2, colony[i]);
    }
}
// 读取城市坐标
void initm()
{
    int i;
    FILE *fp;
    float ffff, eeee;
    if ((fp = fopen(filepath, "r")) == NULL)
    {
        MPI_Finalize();
        exit(-1);
    }
    printf("hello\n");
    fscanf(fp, "%d", &xCity);
    for (i = 0; i < xCity; i++) /*  init cityXY[][]  */
    {
        fscanf(fp, "%*d%f%f", &ffff, &eeee);
        printf("ffff=%f eeee=%f\n", ffff, eeee);
        cityXY[i][0] = ffff;
        cityXY[i][1] = eeee;
    }
    fclose(fp);
}
// 初始化一个种群的染色体
void init()
{
    int i, j, t, sign, mod, array[CITY];
    double d;
    for (i = 0; i < xCity; i++)
    {
        city_dis[i][i] = 0;
        for (j = i + 1; j < xCity; j++)
        {
            d = (cityXY[i][0] - cityXY[j][0]) * (cityXY[i][0] - cityXY[j][0]) * 1.0 +
                (cityXY[i][1] - cityXY[j][1]) * (cityXY[i][1] - cityXY[j][1]) * 1.0;
            city_dis[j][i] = city_dis[i][j] = (int)(sqrt(d) + 0.5);
        }
    }
    mod = xCity;
    for (i = 0; i < xCity; i++)
        array[i] = i; //    init colony[][]
    for (i = 0; i < xColony; i++, mod = xCity)
    {
        for (j = 0; j < xCity; j++)
        {
            sign = rand() % mod;
            colony[i][j] = array[sign];
            t = array[mod - 1];
            array[mod - 1] = array[sign];
            array[sign] = t;
            mod--;
            if (mod == 1)
                colony[i][++j] = array[0];
        }
        //for (j = 0; j < xCity; j++) {
        //  printf("%d ", colony[i][j]);
        //}
        //printf("\n");
    }
    for (i = 0; i < xColony; i++)
    { /*    init dis_p[]       */
        dis_p[i] = 0;
        for (j = 0; j < xCity - 1; j++)
            dis_p[i] += city_dis[colony[i][j]][colony[i][j + 1]];
        dis_p[i] += city_dis[colony[i][0]][colony[i][xCity - 1]];
    }
    ibest = 0;
    sumbest = dis_p[0]; /*  init ibest & sumbest */
    sumTemp = sumbest * 5.0;
    speed = 2147483647;
    loopcounter = 0;
    Ni = 0; /*   initialize GunNum & Ni    */
    //printf("init success!!!\n");
}