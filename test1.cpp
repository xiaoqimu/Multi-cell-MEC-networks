#include<iostream>
using namespace std;
#include<vector>
#include<string>
#include<math.h>
#include <numeric>
#include<set>
#include <cstdlib> 
#include <ctime>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <cstring>
#include <cstdio>
#include <time.h>//clock()需要的头文件
#include <sys/time.h>//gettimeofday()需要的头文件
#include <iomanip>
#include<unordered_map>
#include <cmath>
#include <float.h>
#include <thread> 
#define log2(n) (log(n)/log(2))
vector<vector<vector<vector<double>>>> res111;
vector<vector<vector<double>>> D;
vector<vector<double>> kn;
vector<vector<double>> w_e;
vector<vector<double>> w_t;
vector<vector<double>> w_c;
vector<vector<vector<double>>> h_s;
vector<vector<vector<double>>> I;
vector<vector<pair<double,double>>> f_range_local;
vector<vector<pair<double,double>>> p_range_all;
vector<vector<double>> C;
vector<vector<double>> km;
vector<vector<double>> W_cross;
vector<vector<double>> h_cross;
vector<vector<double>> W;
vector<vector<double>> cost;
vector<vector<double>> I_cross;
double aa;
double cost_u=0;
vector<vector<double>> F_MAX;
vector<vector<pair<double,double>>> f_range_mec;
vector<vector<double>> a;
vector<vector<vector<double>>> a1;
vector<vector<vector<double>>> S;
vector<vector<vector<double>>> V;
vector<vector<double>> a_d;
vector<vector<vector<double>>> a1_d;
vector<vector<vector<double>>> S_d;
vector<vector<vector<double>>> V_d;
vector<vector<vector<double>>> prr;
double noise=pow(10,-14);
double D_MAX=pow(10,10);
long long seed=2;
double percise=pow(10,-10);
double percise1=pow(10,-10);
double fibonacci_itr=pow(10,-20)/6.0;
double v=0.4;
double p=0.1;
double k1=pow(10,-26),k2=pow(10,-25);
int popsize = 1;
int max_iter = 100;
double w1 = 0.6; // 惯性权重
double c1 = 0.5;
double c2 = 1.5;
int Area;
vector<double> fibonacci_num;
vector<double> M;
vector<double> ML;
vector<double> N;
vector<vector<vector<double>>> R1;
vector<vector<vector<double>>> B;
vector<vector<vector<double>>> R1_d;
vector<vector<vector<double>>> B_d;
vector<vector<vector<vector<double>>>> C_V;//0:count;1:value
vector<vector<vector<double>>> S1;//0:count;1:value
vector<vector<vector<double>>> S2;//0:count;1:value
const int MAXN = 10000;
const double eps = 1e-6;
const double INF =pow(10,20);
int lenn; 
int lenm;  
double E_val[MAXN][MAXN];   // 记录每条边的权值
double L_val[MAXN];      // 每个左顶点的值
double R_val[MAXN];       // 每个右顶点的值
bool L_vis[MAXN];    // 记录每一轮匹配匹配过的左顶点
bool R_vis[MAXN];     // 记录每一轮匹配匹配过的右顶点
int match[MAXN];        // 记录每个右顶点匹配到的左顶点 如果没有则为-1
double slack[MAXN];        // 记录每个右顶点如果能被左顶点匹配最少还需要多少值
int FL;
double  f_equal_p_x(int area,int i,int area1,int  j,double p_min,double p_max,double C);
vector<double> f_equal_1(int area,int i,int area1,int j,double x_min,double x_max,int in_de_creasing,int xiao_da,double (*func)(int,int,int,int,double),double C);
vector<double> subpro(int area,int i,int area1,int j);
vector<double> f_equal(int area,int i,int j,double x_min,double x_max,int in_de_creasing,int xiao_da,double (*func)(int,int,int,double),double C);
vector<double> f1_equal(int area,int i,int j,double f_low,double f_high,double x_min,double x_max);
vector<double> f1_equal_1(int area,int i,int area1,int j,double f_low,double f_high,double p_min,double p_max);
template <class T>
void Print(T &x){
    for(auto i:x){
        cout<<i<<' ';
    }
    cout<<endl;
}


bool dfs(int now)
{
    L_vis[now] = true;
    for (int i = 0; i < lenm; i++)
    {
        if (R_vis[i]) continue; // 每一轮匹配 每个右顶点只尝试一次
        double tmp = L_val[now] + R_val[i] - E_val[now][i];
 
        if (fabs(tmp) < eps) {  // 如果符合要求
            R_vis[i] = true;
            if (match[i] == -1 || dfs( match[i] )) {    // 找到一个没有匹配的右顶点 或者该右顶点当前匹配的左顶点可以找到其它匹配
                match[i] = now;
                return true;
            }
        } else {
            slack[i] = min(slack[i], tmp);  // slack 可以理解为该右顶点能被一个左顶点匹配 还需多少值 取最小值
        }
    }
    return false;
}
 
vector<vector<int>> KM(double &opt)
{
    vector<vector<int>> res(lenn,vector<int>(lenm));
    memset(match, -1, sizeof(match));    // 初始每个右顶点都没有匹配的左顶点
    memset(R_val, 0, sizeof(R_val));   // 初始每个右顶点的值为0
    for (int i = 0; i < lenm; ++i) {             //初始每个左顶点为与其相连边的最大权值
        L_val[i] = E_val[i][0];
        for (int j = 1; j < lenm; ++j) {
            L_val[i] = max(L_val[i], E_val[i][j]);
        }
    }
    // 尝试为每一个左顶点匹配
    for (int i = 0; i < lenm; ++i) {
        
        fill(slack, slack + lenm, INF);    // 因为要取最小值 初始化为无穷大
        while (1) {
            // 为每个左顶点匹配的方法是 ：如果找不到就降低期望值，直到找到为止
            // 记录每轮匹配中左右顶点是否被尝试匹配过
            memset(L_vis, false, sizeof(L_vis));
            memset(R_vis, false, sizeof(R_vis));
            if (dfs(i)) break;  // 找到匹配 退出
            // 如果不能找到 就降低期望值
            // 最小可降低的期望值
            double d = INF;
            for (int j = 0; j < lenm; ++j)
                if (!R_vis[j]) d = min(d, slack[j]);
            for (int j = 0; j < lenm; ++j) {
                // 所有访问过的(被涉及的)左顶点降低值
                if (L_vis[j]) L_val[j] -= d;
                // 所有访问过(被涉及的)的右顶点增加值
                if (R_vis[j]) R_val[j] += d;
                // 没有访问过的右顶点 因为左顶点的期望值降低，距离被左顶点匹配又进了一步
                else slack[j] -= d;
            }
        }
    }
    // 匹配完成 求出所有匹配的权值和
    opt = 0;
    int count=0;
    for (int i = 0; i < lenm; ++i)
    if(E_val[match[i]][i]!=0){
        res[match[i]][i]=1;
        count++;
        opt += E_val[match[i]][i];
        //cout<<match[i]<<' '<<i<<' '<<E_val[ match[i] ][i]<<endl;
    }
    //cout<<count<<endl;
    return res;
}


vector<vector<double>> lognorm(double u, double sigma, int row,int column){
    vector<vector<double>> res(row,vector<double>(column));
    mt19937 rand_num(seed);  // 大随机数
	lognormal_distribution<> dist(u, sigma);  // 给定范围	
	for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            res[i][j]= dist(rand_num);   
        }
    }
    seed+=1;
	return res;
}

vector<vector<double>> uniform2(double a, double b, int row,int column){
    vector<vector<double>> res(row,vector<double>(column));
    mt19937 rand_num(seed);  // 大随机数
	uniform_real_distribution< > dist(a, b);  // 给定范围	   
	for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            res[i][j]= dist(rand_num);   
        }
    }
    seed+=1;
	return res;
}
vector<vector<int>> uniform2_int(int a, int b, int row,int column){
    vector<vector<int>> res(row,vector<int>(column));
    mt19937 rand_num(seed);  // 大随机数
	uniform_int_distribution< > dist(a, b);  // 给定范围	   
	for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            res[i][j]= dist(rand_num);   
        }
    }
    seed+=1;
	return res;
}
vector<int> uniform1_int(int a, int b,int column){
    vector<int> res(column);
     mt19937 rand_num(seed);  // 大随机数
	uniform_int_distribution< > dist(a, b);  // 给定范围
    for(int i=0;i<column;i++)
	   res[i]= dist(rand_num);
    seed+=1;
	return res;
}

vector<double> uniform1(double a, double b,int column){
    vector<double> res(column);
     mt19937 rand_num(seed);  // 大随机数
	uniform_real_distribution< > dist(a, b);  // 给定范围
    for(int i=0;i<column;i++)
	   res[i]= dist(rand_num);
    seed+=1;
	return res;
}


vector<double> poisson(int a,int column){
    vector<double> res(column);
     mt19937 rand_num(seed);  // 大随机数
	poisson_distribution<> dist(a);  // 给定范围
    for(int i=0;i<column;i++)
	   res[i]= dist(rand_num);
    seed+=1;
	return res;
}




void fibonacci_creat(){
    fibonacci_num.push_back(1.0);
    fibonacci_num.push_back(1.0);
    double v1=1.0;
    double v2=1.0;
    while(v2<886770.0){
        fibonacci_num.push_back(v1+v2);
        double temp=v1;
        v1=v2;
        v2=temp+v2;
    }

}




void initialization_mec(){
    vector<vector<double>> d(Area);
    ifstream infile("d1.csv");
    string line;
    getline(infile, line);
    for (int i=0;getline(infile, line);i++)
    {
        istringstream sin(line);
        string field;
        getline(sin, field, ',');
        while (getline(sin, field, ',')) 
        {
            d[i].push_back(stod(field.c_str()));
        }
       
    }
    cost.clear();
    M.clear();
    ML.clear();
    f_range_mec.clear();
    W.clear();
    km.clear();
    W_cross.clear();
    h_cross.clear();
    cost.resize(Area);
    M.resize(Area,-1);
    ML.resize(Area,0);
    f_range_mec.resize(Area);
    W.resize(Area);
    km.resize(Area);
    W_cross.resize(Area,vector<double>(Area));
    h_cross.resize(Area,vector<double>(Area));


    for(int i=0;i<Area;i++){
        string name;
        name=name+"d"+to_string(i+2)+".csv";
        ifstream infile(name);
        string line;
        getline(infile, line);
        for (int j=0;getline(infile, line);j++)
        {
            M[i]+=1;
            ML[i]+=1;
            istringstream sin(line);
            string field;
            getline(sin, field, ',');
            getline(sin, field, ',');
            double low_temp=stod(field.c_str())*pow(10,9);
            getline(sin, field, ',');
            double up_temp=stod(field.c_str())*pow(10,9);
            f_range_mec[i].push_back({low_temp,up_temp});
            getline(sin, field, ',');
            double temp_w=stod(field.c_str());
            W[i].push_back(temp_w*pow(10,9));
            getline(sin, field, ',');
            double temp_cost=stod(field.c_str())/3600.0;
            cost[i].push_back(temp_cost);
            
        }
        for(int j=0;j<ML[i];j++){
            km[i].push_back(pow(10,-26));
        }
        
        
    }
    auto temp=lognorm(0.0,pow(8,0.5),Area,Area);


    for(int i=0;i<Area;i++){
        for(int j=0;j<Area;j++){
            if(i!=j){
            
            h_cross[i][j]=128.1+37.6*log10(d[i][j])+temp[i][j];
            W_cross[i][j]=min(W[i][0],W[j][0]);
            }
        }
    }


}

void initialization_user(vector<pair<double,double>> f_local,vector<pair<double,double>> p_all,int lf,int lp){
    kn.clear();
    kn.resize(Area);
    f_range_local.clear();
    f_range_local.resize(Area);
    p_range_all.clear();
    p_range_all.resize(Area);
    D.clear();
    D.resize(Area);
    w_e.clear();
    w_t.clear();
    w_c.clear();
    w_e.resize(Area);
    w_c.resize(Area);
    w_t.resize(Area);
    C.clear();
    C.resize(Area);
    h_s.clear();
    h_s.resize(Area);
    I.clear();
    I.resize(Area);
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        int temp_M=M[i];
        for(int j=0;j<temp_N;j++){
            kn[i].push_back(pow(10,-26));
        }
        auto temp_w=uniform2(0,1,3,temp_N);
        auto temp_D_2=uniform1(1*pow(2,30),10*pow(2,30),temp_N);
        auto temp_D_1=uniform1(1,1.5,temp_N);
        auto temp_D_3=uniform1(0.8,1.2,temp_N);
        auto temp_D_4=uniform1(50,100,temp_N);
        auto temp_D_5=uniform1(10800,18000,temp_N);
        auto temp_D_6=uniform1(pow(10,7),pow(10,8),temp_N);
        auto temp_C=uniform1(1,2,temp_N);
        auto temp_h=lognorm(0.0,pow(8,0.5),temp_N,temp_M+1);
        auto temp_d1=uniform1(1,100,temp_N);
        auto temp_d2=uniform2(1,10,temp_N,temp_M);
        auto temp_I=uniform2(0,30,temp_N,temp_M);
        for(int j=0;j<temp_N;j++){
            double sum_temp_w=temp_w[0][j]+temp_w[1][j]+temp_w[2][j];
            temp_w[0][j]=temp_w[0][j]/sum_temp_w;
            temp_w[1][j]=temp_w[1][j]/sum_temp_w;
            temp_w[2][j]=temp_w[2][j]/sum_temp_w;
            w_t[i].push_back(temp_w[0][j]);
            w_c[i].push_back(temp_w[1][j]/cost_u);
            w_e[i].push_back(temp_w[2][j]/aa);
           
            srand(seed);
            int index_f=rand()%(lf)+0;
            
            int index_p=rand()%(lp)+0;
            seed+=1;
            f_range_local[i].push_back(f_local[index_f]);
            p_range_all[i].push_back(p_all[index_p]);
            
            D[i].push_back({temp_D_1[j]*temp_D_2[j],temp_D_2[j],temp_D_3[j]*temp_D_2[j],temp_D_4[j],temp_D_5[j],temp_D_6[j]});
            C[i].push_back(temp_C[j]*aa*(18000+pow(10,8)/aa)/(temp_D_5[j]+temp_D_6[j]/aa));
            vector<double> h_s_temp;
            vector<double> I_temp;
            I_temp.push_back(0);
            h_s_temp.push_back(128.1 + 37.6*log10(temp_d1[j]) + temp_h[j][0]);
            for(int jj=1;jj<temp_M+1;jj++){
                h_s_temp.push_back(140.7 + 36.7*log10(temp_d2[j][jj-1]) + temp_h[j][jj]);
                I_temp.push_back(temp_I[j][jj-1]);
            }
            h_s[i].push_back(h_s_temp);
            I[i].push_back(I_temp);
        }

    }

}



double func_local(int area,int i,double f){
    return w_t[area][i]*D[area][i][1]*D[area][i][3]/f+kn[area][i]*w_e[area][i]*D[area][i][1]*D[area][i][3]*f*f;
}

double func_MBS_SBS(int area,int i, int j,double x,double f){
   double res= (D[area][i][0] + D[area][i][2]) * ((w_t[area][i]+w_c[area][i]*cost[area][j]) * x + w_e[area][i] * (noise+I[area][i][j]) * x * (pow(2.0,1.0 / (W[area][j] * x)) - 1) / h_s[area][i][j]) +(w_t[area][i]+w_c[area][i]*cost[area][j]) * D[area][i][1] * D[area][i][3] / f + km[area][j] * w_e[area][i] * D[area][i][1] * D[area][i][3] * f * f;
   return res;
}


double func_R_1(int area,int i,int area1,int j,double p){
    return 1 / (W[area][0] * (log2(1 + p * h_s[area][i][0] / (noise+I[area][i][0]))))+1 / (W_cross[area][area1] * (log2(1 + p * h_cross[area][area1] / (noise))));
}



double func_MBS_SBS_1(int area,int i,int area1,int j,double p,double f){
    return (D[area][i][0] + D[area][i][2]) * (w_t[area][i]+w_c[area][i]*cost[area1][j]+w_e[area][i]*p)*func_R_1(area,i,area1,j,p)+(w_t[area][i]+w_c[area][i]*cost[area1][j]) * D[area][i][1] * D[area][i][3] / f + km[area1][j] * w_e[area][i] * D[area][i][1] * D[area][i][3] * f * f;
}



 
double func_df(int area,int i, int j,double x){
    return (w_t[area][i]+w_c[area][i]*cost[area][j]) + w_e[area][i] * (noise+I[area][i][j]) * (pow(2.0,1.0 / (W[area][j] * x)) * (1 - log(2) / (x * W[area][j])) - 1) / h_s[area][i][j];
}




double func_df_1(int area,int i,int area1,int j,double p){
    return  w_t[area][i]+w_c[area][i]*cost[area1][j] + w_e[area][i]*p+w_e[area][i]*func_R_1(area,i,area1,j,p)/(-h_s[area][i][0]/((noise+I[area][i][0]+h_s[area][i][0]*p)*W[area][0]*pow(log2(1+h_s[area][i][0]*p/(noise+I[area][i][0])),2.0))-h_cross[area][area1]/((noise+h_cross[area][area1]*p)*W_cross[area][area1]*pow(log2(1+h_cross[area][area1]*p/(noise)),2.0)));
}




double func_f1(int area,int i,int j,double x){
    return (noise+I[area][i][j]) * x * (pow(2.0,1.0 / (W[area][j] * x)) - 1) / h_s[area][i][j];
}



double func_f1_p(int area,int i,int area1,int j,double p){
    return p*(1 / (W[area][0] * (log2(1 + p * h_s[area][i][0] / (noise+I[area][i][0]))))+1 / (W_cross[area][area1] * (log2(1 + p * h_cross[area][area1] / (noise)))));
}


double func_R(int area,int i,int  j,double p){
    return 1 / (W[area][j] * (log2(1 + p * h_s[area][i][j] / (noise+I[area][i][j]))));
}
double func_MBS_SBS_11(int area,int i,int area1,int j,double p,double f){
    if(area==area1){
        return (D[area][i][0] + D[area][i][2]) * (w_t[area][i]+w_c[area][i]*cost[area][j]+w_e[area][i]*p)*func_R(area,i,j,p)+(w_t[area][i]+w_c[area][i]*cost[area][j]) * D[area][i][1] * D[area][i][3] / f + km[area][j] * w_e[area][i] * D[area][i][1] * D[area][i][3] * f * f;
    }
    else{
        return (D[area][i][0] + D[area][i][2]) * (w_t[area][i]+w_c[area][i]*cost[area1][j]+w_e[area][i]*p)*func_R_1(area,i,area1,j,p)+(w_t[area][i]+w_c[area][i]*cost[area1][j]) * D[area][i][1] * D[area][i][3] / f + km[area1][j] * w_e[area][i] * D[area][i][1] * D[area][i][3] * f * f;
    }
    
}

double func_MBS_SBS_111(int area,int i,int area1,int j,double p,double f,vector<vector<double>> x){
    if(area==area1){
        return (D[area][i][0] + D[area][i][2]) * (w_t[area][i]+w_c[area][i]*cost[area][j]+w_e[area][i]*p)*func_R(area,i,j,p)+(w_t[area][i]+w_c[area][i]*cost[area][j]) * D[area][i][1] * D[area][i][3] / f + km[area][j] * w_e[area][i] * D[area][i][1] * D[area][i][3] * f * f;
    }
    else{
        return (D[area][i][0] + D[area][i][2]) * (w_t[area][i]+w_c[area][i]*cost[area1][j]+w_e[area][i]*p)*(func_R(area,i,0,p)+x[area][area1])+(w_t[area][i]+w_c[area][i]*cost[area1][j]) * D[area][i][1] * D[area][i][3] / f + km[area1][j] * w_e[area][i] * D[area][i][1] * D[area][i][3] * f * f;
    }
    
}

double func_R_rev(int area,int i,int j,double x){//单元内
    return (noise+I[area][i][j]) * (pow(2.0,1.0 / (W[area][j] * x)) - 1) / h_s[area][i][j];
}


double func_df1_mix(int area,int i,int  j,double x){
    return (D[area][i][0] + D[area][i][2]) * (noise+I[area][i][j]) * (pow(2.0,1.0 / (W[area][j] * x)) * (1 - log(2) / (x * W[area][j])) - 1) / h_s[area][i][j]+ (2 * (D[area][i][0] + D[area][i][2]) * km[area][j] * pow(D[area][i][1] * D[area][i][3],3.0)) / (pow(D[area][i][4] - (D[area][i][0]+ D[area][i][2]) * x,3.0));
}



double func_f1_mix(int area,int i,int j,double x){
    return (D[area][i][0] + D[area][i][2]) * func_f1(area,i, j,x) + (km[area][j] * pow(D[area][i][1] * D[area][i][3],3.0))/(pow(D[area][i][4] - (D[area][i][0] + D[area][i][2]) * x,2.0)) - D[area][i][5];
}



double func_df1_mix_1(int area,int i,int area1,int j,double p){
    return (D[area][i][0] + D[area][i][2]) *(func_R_1(area,i,area1,j,p)/(-h_s[area][i][0]/((noise+I[area][i][0]+h_s[area][i][0]*p)*W[area][0]*pow(log2(1+h_s[area][i][0]*p/(noise+I[area][i][0])),2.0))-h_cross[area][area1]/((noise+h_cross[area][area1]*p)*W_cross[area][area1]*pow(log2(1+h_cross[area][area1]*p/(noise)),2.0)))+p)+ (2 * (D[area][i][0] + D[area][i][2]) * km[area1][j] * pow(D[area][i][1] * D[area][i][3],3.0)) / (pow(D[area][i][4] - (D[area][i][0]+ D[area][i][2]) * func_R_1(area,i,area1,j,p),3.0));
}



double func_f1_mix_1(int area,int i,int area1,int  j,double p){
    return (D[area][i][0] + D[area][i][2]) *func_R_1(area,i,area1,j,p)*p+ (km[area1][j] * pow(D[area][i][1] * D[area][i][3],3.0)) /(pow(D[area][i][4] - (D[area][i][0] + D[area][i][2]) * func_R_1(area,i,area1,j,p),2.0)) - D[area][i][5];
}


vector<double> f_to_x_opt_MBS(int area,int i){
    vector<double> res{0,0};
    double x_MBS_min = func_R(area,i,0, p_range_all[area][i].second);
    double x_MBS_max = func_R(area,i,0, p_range_all[area][i].first);
    //cout<<x_MBS_min<<' '<<x_MBS_max<<endl;
    auto temp_res1 = f1_equal(area,i,0,f_range_mec[area][0].first,f_range_mec[area][0].first,x_MBS_min,x_MBS_max);
    //cout<<temp_res1[0]<<' '<<temp_res1[1]<<' '<<temp_res1[2]<<endl;
    res[0]=temp_res1[0];    
    if(res[0]!=0){
        if(func_df(area,i,0,temp_res1[1])>=0){
            res[1]=temp_res1[1];
        }
        else if(func_df(area,i,0,temp_res1[2])<=0){
            res[1]=temp_res1[2];
        }
        else{
            double low = temp_res1[1],high = temp_res1[2],mid = 0.0;
            while(low <= high){
                mid=(low+high)/2.0;
                double temp=func_df(area,i,0,mid);
                if(abs(temp)<=percise){
                    break;
                }
                else{
                    if(temp - percise >= 0){
                        high=mid-percise1;
                    }
                    else{
                        low=mid+percise1;
                    }
                }
            }
            res[1]=mid;
        }
    }
    return res;
}



vector<double> f_to_x_opt_MBS_1(int area,int i,int area1){
    vector<double> res{0,0};
    double p_MBS_min =  p_range_all[area][i].first;
    double p_MBS_max = p_range_all[area][i].second;
    
    auto temp_res1 = f1_equal_1(area,i,area1,0, f_range_mec[area][0].first, f_range_mec[area][0].first, p_MBS_min,p_MBS_max);
    res[0]=temp_res1[0];
    if(temp_res1[0]!=0){
        if(func_df_1(area,i,area1,0,temp_res1[2])>=0){
            res[1]=temp_res1[2]; 
        }
            
       else if(func_df_1(area,i,area1,0,temp_res1[1])<=0){
            res[1]=temp_res1[1];
       }
       else{
            double low = temp_res1[1],high = temp_res1[2],mid = 0.0;
            while(low <= high){
                mid = (low + high) / 2.0;
                double temp=func_df_1(area,i,area1,0,mid);
                if(abs(temp) <= percise){
                    break;
                }
                else{
                    if(temp - percise > 0){
                        low=mid+percise1;
                    }
                        
                    else{
                        high=mid-percise1;
                    }
                       
                }
            }
            res[1]=mid;
       }
            
    }
    return res;
}

vector<double> f_to_x_opt_MBS1(int area,int i,int j,double f){//本单元边缘服务器上，将f视为常数求解p
    vector<double> res{0,0};
    double x_MBS_min = func_R(area,i,j, p_range_all[area][i].second);
    double x_MBS_max = func_R(area,i,j, p_range_all[area][i].first);
    //cout<<x_MBS_min<<' '<<x_MBS_max<<endl;
    auto temp_res1 = f1_equal(area,i,j,f,f,x_MBS_min,x_MBS_max);
    //cout<<temp_res1[0]<<' '<<temp_res1[1]<<' '<<temp_res1[2]<<endl;
    res[0]=temp_res1[0];    
    if(res[0]!=0){
        if(func_df(area,i,j,temp_res1[1])>=0){
            res[1]=temp_res1[1];
        }
        else if(func_df(area,i,j,temp_res1[2])<=0){
            res[1]=temp_res1[2];
        }
        else{
            double low = temp_res1[1],high = temp_res1[2],mid = 0.0;
            while(low <= high){
                mid=(low+high)/2.0;
                double temp=func_df(area,i,j,mid);
                if(abs(temp)<=percise){
                    break;
                }
                else{
                    if(temp - percise >= 0){
                        high=mid-percise1;
                    }
                    else{
                        low=mid+percise1;
                    }
                }
            }
            res[1]=mid;
        }
    }
    return res;
}



vector<double> f_to_x_opt_MBS1_1(int area,int i,int area1,int j,double f){
    vector<double> res{0,0};
    double p_MBS_min =  p_range_all[area][i].first;
    double p_MBS_max = p_range_all[area][i].second;
    auto temp_res1 = f1_equal_1(area,i,area1,j, f, f, p_MBS_min,p_MBS_max);
    res[0]=temp_res1[0];
    if(temp_res1[0]!=0){
        if(func_df_1(area,i,area1,j,temp_res1[2])>=0){
            res[1]=temp_res1[2]; 
        }
            
       else if(func_df_1(area,i,area1,j,temp_res1[1])<=0){
            res[1]=temp_res1[1];
       }
       else{
            double low = temp_res1[1],high = temp_res1[2],mid = 0.0;
            while(low <= high){
                mid = (low + high) / 2.0;
                double temp=func_df_1(area,i,area1,j,mid);
                if(abs(temp) <= percise){
                    break;
                }
                else{
                    if(temp - percise > 0){
                        low=mid+percise1;
                    }
                        
                    else{
                        high=mid-percise1;
                    }
                       
                }
            }
            res[1]=mid;
       }
            
    }
    return res;
}


vector<double> f_to_x_range(int area,int i,int j,double x_min,double x_max){
    vector<double> res{0,0,0};
    res = f1_equal(area,i,j,f_range_mec[area][j].first, f_range_mec[area][j].second, x_min, x_max);
    if(res[0]!=0){
        if(func_df1_mix(area,i,j,res[1]) >= 0){
            res = f_equal(area,i,j,res[1],res[2], 1, 1, func_f1_mix, 0);
            if(res[0]==0){
                return res;
            }
        }
        else if(func_df1_mix(area,i, j,res[2]) <= 0){
            res=f_equal(area,i, j,res[1], res[2], 0, 1, func_f1_mix, 0);
            if(res[0]==0){
                return res;
            }
        }
        else{
            double low = res[1],high = res[2],mid = 0.0;
            while(low <= high){
                mid = (low + high) / 2.0;
                auto temp=func_df1_mix(area,i,j,mid);
                if(abs(temp) <= percise)
                    break;
                else{
                    if(temp- percise > 0)
                        high=mid-percise1;
                    else
                        low=mid+percise1;
                }
            }
            auto temp_res1= f_equal(area,i, j,res[1],mid, 0, 1,func_f1_mix, 0);
            res[0]=temp_res1[0];
            res[1]=temp_res1[1];
            if(res[0]==0){
                return res;
            }
            auto temp_res2= f_equal(area,i,j,mid,res[2], 1, 1,func_f1_mix, 0);
            res[0]=temp_res2[0];
            res[2]=temp_res2[2];
            if(res[0]==0){
                return res;
            }
        }
    }

    return res;

}


vector<double> f_to_x_range_1(int area,int i,int area1,int j,double x_min,double x_max){
    vector<double> res{0,0,0};
    res = f1_equal_1(area,i,area1,j,f_range_mec[area1][j].first, f_range_mec[area1][j].second, x_min, x_max);
    if(res[0]!=0){
        if(func_df1_mix_1(area,i,area1,j,res[2]) >= 0){
            res = f_equal_1(area,i,area1,j,res[1],res[2], 0, 1,func_f1_mix_1, 0);
            if(res[0]==0){
                return res;
            }
        }
        else if(func_df1_mix_1(area,i,area1, j,res[1]) <= 0){
            res=f_equal_1(area,i,area1, j,res[1], res[2], 1, 1, func_f1_mix_1, 0);
            if(res[0]==0){
                return res;
            }
        }
        else{
            double low = res[1],high = res[2],mid = 0.0;
            while(low <= high){
                mid = (low + high) / 2.0;
                auto temp=func_df1_mix_1(area,i,area1,j,mid);
                if(abs(temp) <= percise)
                    break;
                else{
                    if(temp- percise > 0)
                        low=mid+percise1;
                    else
                        high=mid-percise1;
                }
            }
            auto temp_res1= f_equal_1(area,i, area1,j,res[1],mid, 0, 1,func_f1_mix_1, 0);
            res[0]=temp_res1[0];
            res[1]=temp_res1[1];
            if(res[0]==0){
                return res;
            }
            auto temp_res2= f_equal_1(area,i,area1,j,mid,res[2], 1, 1,func_f1_mix_1, 0);
            res[0]=temp_res2[0];
            res[2]=temp_res2[2];
            if(res[0]==0){
                return res;
            }
        }
    }

    return res;
}


double x_to_f_opt(int area,int i,int j,double x){
    double res=0;
    double f_min = max(f_range_mec[area][j].first,(D[area][i][1]*D[area][i][3]) / (D[area][i][4] - (D[area][i][0] + D[area][i][2]) * x));
    double f_max = min(f_range_mec[area][j].second,sqrt((D[area][i][5] - ((D[area][i][0] + D[area][i][2]) * func_f1(area,i, j,x))) / (km[area][j] * D[area][i][1] * D[area][i][3])));
    if(f_min<=f_max){
        double f_opt=pow((w_t[area][i]+w_c[area][i]*cost[area][j])/(2 * w_e[area][i] * km[area][j]),1.0/3.0);
        if(f_opt<f_min){
            res=f_min;
        }
        else if(f_opt>f_max){
            res=f_max;
        }
        else{
            res=f_opt;
        }
    }
    return res;
}

double x_to_f_opt_1(int area,int i,int area1,int j,double x,double p){
    double res=0;
    double f_min = max(f_range_mec[area1][j].first,(D[area][i][1]*D[area][i][3]) / (D[area][i][4] - (D[area][i][0] + D[area][i][2]) * x));
    double f_max = min(f_range_mec[area1][j].second, sqrt((D[area][i][5] - ((D[area][i][0] + D[area][i][2])*x*p)) / (km[area1][j] * D[area][i][1] * D[area][i][3])));
    if(f_min<=f_max){
            double f_opt=pow((w_t[area][i]+w_c[area][i]*cost[area1][j])/(2 * w_e[area][i] * km[area1][j]),1.0/3.0);
            if(f_opt<f_min){
                res=f_min;
            }
            else if(f_opt>f_max){
                res=f_max;
            }
            else{
                res=f_opt;
            }
    }
    return res;
}




vector<double> fibonacci_search(int area,int i,int j,double x_min,double x_max,double c){
    vector<double> res{0,0,0};//f,x,g
    int N=FL;
    int k=N-3;
    double w = 0.0,a,b,f_opt_a,f_opt_b,f_a,f_b;
    a = x_min + (x_max - x_min) * (fibonacci_num[N-1] - fibonacci_num[N-2]) / fibonacci_num[N-1];
    b = x_min + (x_max - x_min) * fibonacci_num[N-2] / fibonacci_num[N-1];
    while(k>=1){
        f_opt_a = x_to_f_opt(area,i, j,a);
        f_a = func_MBS_SBS(area,i, j,a, f_opt_a); 
        f_opt_b = x_to_f_opt(area,i, j,b); 
        f_b = func_MBS_SBS(area,i, j,b, f_opt_b);
        if(k==1){
            w=c;
        }
        if(f_a<f_b){
            x_max=b;
            b=a;
            a = x_min + (x_max - x_min) * ((fibonacci_num[k+1] - fibonacci_num[k]) / fibonacci_num[k+1] - w);
        }
            
        else{
            x_min = a;
            a = b;
            b = x_min + (x_max - x_min) * (fibonacci_num[k] / fibonacci_num[k+1] + w);
        }
        k=k-1;
    }
    f_opt_a = x_to_f_opt(area,i, j,a);
    f_a = func_MBS_SBS(area,i, j,a, f_opt_a);
    f_opt_b = x_to_f_opt(area,i, j,b);
    f_b = func_MBS_SBS(area,i, j,b, f_opt_b);
    if(f_a < f_b){
        x_max = b;
    }
    else{
        x_min = a;
    }
    res[1] = (x_max + x_min) / 2.0;   
    res[0] = x_to_f_opt(area,i, j,res[1]);
    res[2] = func_MBS_SBS(area,i, j,res[1], res[0]);
    return res;
}


vector<double> fibonacci_search_1(int area,int i,int area1,int j,double p_min,double  p_max,double c){
    vector<double> res{0,0,0};//f,p,g
    int N = FL,k=N-3;
    double  w = 0.0,x_min,x_max,a,b,pa,pb,f_opt_a,f_a,f_b,f_opt_b;

    x_min=func_R_1(area,i,area1,j,p_max);
    x_max=func_R_1(area,i,area1,j,p_min);
    a = x_min + (x_max - x_min) * (fibonacci_num[N-1] - fibonacci_num[N-2]) / fibonacci_num[N-1];
    b = x_min + (x_max - x_min) * fibonacci_num[N-2] / fibonacci_num[N-1];
    
    while(k>=1){
        pa= f_equal_p_x(area,i,area1, j,p_min,p_max, a);
        pb=f_equal_p_x(area,i,area1, j,p_min,p_max, b);
        f_opt_a = x_to_f_opt_1(area,i,area1, j,a,pa);
        f_a = func_MBS_SBS_1(area,i,area1, j,pa, f_opt_a);
        f_opt_b = x_to_f_opt_1(area,i,area1, j,b,pb);
        f_b = func_MBS_SBS_1(area,i,area1, j,pb, f_opt_b);
        if(k==1){
            w=c;
        }
        if(f_a<f_b){
            x_max=b;
            b=a;
            a = x_min + (x_max - x_min) * ((fibonacci_num[k+1] - fibonacci_num[k]) / fibonacci_num[k+1] - w);
        }
        else{
            x_min = a;
            a = b;
            b = x_min + (x_max - x_min) * (fibonacci_num[k] / fibonacci_num[k+1] + w);
        }
        k=k-1;
    }
    pa= f_equal_p_x(area,i,area1, j,p_min,p_max, a);
    pb=f_equal_p_x(area,i,area1, j,p_min,p_max, b);
    f_opt_a = x_to_f_opt_1(area,i,area1, j, a,pa);
    f_a = func_MBS_SBS_1(area,i, area1,j, pa, f_opt_a);
    f_opt_b = x_to_f_opt_1(area,i,area1, j,b,pb);
    f_b = func_MBS_SBS_1(area,i,area1, j,pb, f_opt_b);
    if(f_a < f_b){
        x_max = b;
    }
    else{
        x_min = a;
    }   
    double x_opt = (x_max + x_min) / 2.0;
   // cout<<x_opt<<endl;
    res[1]= f_equal_p_x(area,i,area1, j,p_min,p_max, x_opt);
   // cout<<res[1]<<endl;
    res[0] = x_to_f_opt_1(area,i,area1, j, x_opt,res[1]);
    res[2] = func_MBS_SBS_1(area,i, area1,j, res[1], res[0]);
    //Print(res);
    return res;
}

double  f_equal_p_x(int area,int i,int area1,int  j,double p_min,double p_max,double C){
    double low = p_min, high = p_max, mid = 0.0;
    while(low<=high){
        mid = (low + high) / 2.0;
        double temp=func_R_1(area,i,area1,j,mid) - C;
        if(abs(temp) <= percise){
            break;
        }
        else{
            if(temp-percise>0)
                low=mid+percise1;
            else
                high=mid-percise1;
        }
    }
    return mid;
}


vector<double> f1_equal(int area,int i,int j,double f_low,double f_high,double x_min,double x_max){
    vector<double> res{0,0,0};
    double x_max_temp = min(x_max, (D[area][i][4] - D[area][i][1]*D[area][i][3] / f_high) / (D[area][i][0] + D[area][i][2]));
    if(x_min > x_max_temp)
        return res;
    double C = (D[area][i][5] - km[area][j] * D[area][i][1] * D[area][i][3] * f_low * f_low) / (D[area][i][0] + D[area][i][2]);
    res = f_equal(area,i, j,x_min, x_max_temp, 0, 1, func_f1, C);
    return res;
}
vector<double> f1_equal_1(int area,int i,int area1,int j,double f_low,double f_high,double p_min,double p_max){
    vector<double> res{0,0,0};
    double C;
    C=(D[area][i][4] - D[area][i][1]*D[area][i][3] / f_high) / (D[area][i][0] + D[area][i][2]);
    auto temp_res =f_equal_1(area,i,area1, j,p_min, p_max,0,1,func_R_1 , C);
    if(temp_res[0]==0)
        return res;
    C =  (D[area][i][5] - km[area1][j] * D[area][i][1] * D[area][i][3] * f_low * f_low) / (D[area][i][0] + D[area][i][2]);
    res = f_equal_1(area,i,area1,j,temp_res[1],temp_res[2], 1, 1, func_f1_p, C);
    return res;
}

vector<double> f_equal(int area,int i,int j,double x_min,double x_max,int in_de_creasing,int xiao_da,double (*func)(int,int,int,double),double C){
    //in_de_creasing is 1 denote increasing; otherwise
    //xiao_da is 1 denote less than or equal; otherwise
    vector<double> res{0,x_min,x_max};
    
    if(in_de_creasing==1){
        double temp1=(*func)(area,i,j,x_min),temp2=(*func)(area,i,j,x_max);
        if(xiao_da==1){
            if(temp1>C){
                return res;
            }
            else if(temp2<C){
                res[0]=1;
                return res;
            }
        }
        else{
            if(temp2<C){
                return res;
            }
            else if(temp1>C){
                res[0]=1;
                return res;
            }    
        }
    }
    else{
        double temp1=(*func)(area,i,j,x_min),temp2=(*func)(area,i,j,x_max);
        if(xiao_da==1){
            if(temp2>C){
                return res;
            }
            else if(temp1<C){
                res[0]=1;
                return res;
            }
    
        }
        else{
            if(temp1<C){
                return res;
            }
            else if(temp2>C){
                res[0]=1;
                return res;
            }

        }
        
    
    }
    double low = x_min, high = x_max, mid = 0.0;
    while(low<=high){
        mid = (low + high) / 2.0;
        double temp=(*func)(area,i,j,mid)-C;
        if(abs(temp) <= percise)
            break;
        else{
            if(temp-percise>0){
                if(in_de_creasing==1)
                    high=mid-percise1;
                else
                    low=mid+percise1;
            }
            else{
                if(in_de_creasing==1)
                    low=mid+percise1;
                else
                    high=mid-percise1;
            }
        }
        
    }

    if((xiao_da==1&&in_de_creasing==1) || (in_de_creasing==0 && xiao_da==0)){
        res[2] = mid;
    }
    else{
        res[1] = mid;
    }
    res[0]=1;
    return res;
}

vector<double> f_equal_1(int area,int i,int area1,int j,double x_min,double x_max,int in_de_creasing,int xiao_da,double (*func)(int,int,int,int,double),double C){
    vector<double> res{0,x_min,x_max};
    if(in_de_creasing==1){
        double temp1=(*func)(area,i,area1,j,x_min),temp2=(*func)(area,i,area1,j,x_max);
        if(xiao_da==1){
            if(temp1>C){
                return res;
            }
            else if(temp2<C){
                res[0]=1;
                return res;
            }
        }
        else{
            if(temp2<C){
                return res;
            }
            else if(temp1>C){
                res[0]=1;
                return res;
            }    
        }
    }
    else{
        double temp1=(*func)(area,i,area1,j,x_min),temp2=(*func)(area,i,area1,j,x_max);
        if(xiao_da==1){
            if(temp2>C){
                return res;
            }
            else if(temp1<C){
                res[0]=1;
                return res;
            }
    
        }
        else{
            if(temp1<C){
                return res;
            }
            else if(temp2>C){
                res[0]=1;
                return res;
            }

        }
        
    
    }
    double low = x_min, high = x_max, mid = 0.0;
    while(low<=high){
        mid = (low + high) / 2.0;
        double temp=(*func)(area,i,area1,j,mid)-C;
        if(abs(temp) <= percise)
            break;
        else{
            if(temp-percise>0){
                if(in_de_creasing==1)
                    high=mid-percise1;
                else
                    low=mid+percise1;
            }
            else{
                if(in_de_creasing==1)
                    low=mid+percise1;
                else
                    high=mid-percise1;
            }
        }
        
    }

    if((xiao_da==1&&in_de_creasing==1) || (in_de_creasing==0 && xiao_da==0)){
        res[2] = mid;
    }
    else{
        res[1] = mid;
    }
    res[0]=1;
    return res;
}


void subpro1(int area,int i){
    //一个用户需要完成的所有计算
    int M_s=accumulate(M.begin(),M.end(),0);
    vector<double> res_g(M_s+Area+1);//G
    vector<double> res_f(M_s+1);//f
    vector<double> res_p(M_s+Area);//p
    
    double f_local_min=max(f_range_local[area][i].first,D[area][i][1]*D[area][i][3]/D[area][i][4]);
    double f_local_max=min(f_range_local[area][i].second,sqrt(D[area][i][5]/(kn[area][i]*D[area][i][1]*D[area][i][3])));
    if(f_local_min<=f_local_max){
        double f_local_opt=pow(w_t[area][i]/(2*w_e[area][i]*kn[area][i]),1.0/3.0);
        if(f_local_opt<f_local_min){
            res_g[0]=func_local(area,i,f_local_min);
            res_f[0]=f_local_min;
        }
        else if(f_local_opt>f_local_max){
            res_g[0]=func_local(area,i,f_local_max);
            res_f[0]=f_local_max;
        }
        else{
            res_g[0]=func_local(area,i,f_local_opt);
            res_f[0]=f_local_opt;
        }
    }
    int index_g=1,index_p=0,index_f=1;  
    for(int ii=0;ii<Area;ii++){
        int temp_M=M[ii];
        auto subres=subpro(area,i,ii,0);
        res_g[index_g++]=subres[0];
        res_p[index_p++]=subres[2];
        for(int k=0;k<temp_M;k++){
            if(F_MAX[ii][k]==0){
                index_g++;
                index_f++;
                index_p++;
            }
            else{
                auto subres=subpro(area,i,ii,k+1);
                res_g[index_g++]=subres[0];
                res_f[index_f++]=subres[1];
                res_p[index_p++]=subres[2];
            }
                    
        }
    }
    res111[0][area][i]=res_g;
    res111[1][area][i]=res_f;
    res111[2][area][i]=res_p;
}

void solution_11(){
    //一个用户一个线程
    res111.clear();
    res111.resize(3,vector<vector<vector<double>>>(Area));
    int M_s=accumulate(M.begin(),M.end(),0.0);
    vector<thread> thread_vec; 
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        res111[0][i].resize(temp_N,vector<double>(M_s+Area+1));//G
        res111[1][i].resize(temp_N,vector<double>(M_s+1));//f
        res111[2][i].resize(temp_N,vector<double>(M_s+Area));//p
    }
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        for(int j=0;j<temp_N;j++){
            thread_vec.emplace_back(thread(subpro1,i,j));
        }
    }
    for (thread& t : thread_vec) 
       t.join(); 

}

void subpro11(int area){
    //一个用户需要完成的所有计算
    int M_s=accumulate(M.begin(),M.end(),0);
    vector<vector<double>> res_g(N[area],vector<double>(M_s+Area+1));//G
    vector<vector<double>> res_f(N[area],vector<double>(M_s+1));//f
    vector<vector<double>> res_p(N[area],vector<double>(M_s+Area));//p
    for(int i=0;i<N[area];i++){
    double f_local_min=max(f_range_local[area][i].first,D[area][i][1]*D[area][i][3]/D[area][i][4]);
    double f_local_max=min(f_range_local[area][i].second,sqrt(D[area][i][5]/(kn[area][i]*D[area][i][1]*D[area][i][3])));
    if(f_local_min<=f_local_max){
        double f_local_opt=pow(w_t[area][i]/(2*w_e[area][i]*kn[area][i]),1.0/3.0);
        if(f_local_opt<f_local_min){
            res_g[i][0]=func_local(area,i,f_local_min);
            res_f[i][0]=f_local_min;
        }
        else if(f_local_opt>f_local_max){
            res_g[i][0]=func_local(area,i,f_local_max);
            res_f[i][0]=f_local_max;
        }
        else{
            res_g[i][0]=func_local(area,i,f_local_opt);
            res_f[i][0]=f_local_opt;
        }
    }
    int index_g=1,index_p=0,index_f=1;  
    for(int ii=0;ii<Area;ii++){
        int temp_M=M[ii];
        auto subres=subpro(area,i,ii,0);
        res_g[i][index_g++]=subres[0];
        res_p[i][index_p++]=subres[2];
        for(int k=0;k<temp_M;k++){
            if(F_MAX[ii][k]==0){
                index_g++;
                index_f++;
                index_p++;
            }
            else{
                auto subres=subpro(area,i,ii,k+1);
                res_g[i][index_g++]=subres[0];
                res_f[i][index_f++]=subres[1];
                res_p[i][index_p++]=subres[2];
            }
                    
        }
    }
    }
    res111[0][area]=res_g;
    res111[1][area]=res_f;
    res111[2][area]=res_p;
}

void solution_1111(){
    //一个单元一个线程
    res111.clear();
    res111.resize(3,vector<vector<vector<double>>>(Area));
    int M_s=accumulate(M.begin(),M.end(),0.0);
    vector<thread> thread_vec; 
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        res111[0][i].resize(temp_N,vector<double>(M_s+Area+1));//G
        res111[1][i].resize(temp_N,vector<double>(M_s+1));//f
        res111[2][i].resize(temp_N,vector<double>(M_s+Area));//p
    }
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        thread_vec.emplace_back(thread(subpro11,i));
    }
    for (thread& t : thread_vec) 
       t.join(); 

}

 vector<vector<vector<vector<double>>>> solution_1(){
    vector<vector<vector<vector<double>>>> res(3,vector<vector<vector<double>>>(Area));
    int M_s=accumulate(M.begin(),M.end(),0.0);
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        res[0][i].resize(temp_N,vector<double>(M_s+Area+1));//G
        res[1][i].resize(temp_N,vector<double>(M_s+1));//f
        res[2][i].resize(temp_N,vector<double>(M_s+Area));//p
        for(int j=0;j<temp_N;j++){
            int index_g=0,index_p=0,index_f=0;  
            auto subres=subpro(i,j,i,-1);
            res[0][i][j][index_g++]=subres[0];
            res[1][i][j][index_f++]=subres[1];
            for(int ii=0;ii<Area;ii++){
                int temp_M=M[ii];
                auto subres=subpro(i,j,ii,0);
                res[0][i][j][index_g++]=subres[0];
                res[2][i][j][index_p++]=subres[2];
                for(int k=0;k<temp_M;k++){
                    if(F_MAX[ii][k]==0){
                        index_g++;
                        index_f++;
                        index_p++;
                    }
                    else{
                        auto subres=subpro(i,j,ii,k+1);
                        res[0][i][j][index_g++]=subres[0];
                        res[1][i][j][index_f++]=subres[1];
                        res[2][i][j][index_p++]=subres[2];
                    }
                    
                }
            }
        }
    }

    return res;

}


vector<vector<vector<vector<double>>>> solution_111(){
    vector<vector<vector<vector<double>>>> res(3,vector<vector<vector<double>>>(Area));
    int M_s=accumulate(M.begin(),M.end(),0.0);
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        res[0][i].resize(temp_N,vector<double>(M_s+Area+1));//G
        res[1][i].resize(temp_N,vector<double>(M_s+1));//f
        res[2][i].resize(temp_N,vector<double>(M_s+Area));//p
        for(int j=0;j<temp_N;j++){
            int index_g=0,index_p=0,index_f=0;  
            auto subres=subpro(i,j,i,-1);
            res[0][i][j][index_g++]=subres[0];
            res[1][i][j][index_f++]=subres[1];
            for(int ii=0;ii<Area;ii++){
                int temp_M=M[ii];
                auto subres=subpro(i,j,ii,0);
                res[0][i][j][index_g++]=subres[0];
                res[2][i][j][index_p++]=subres[2];
                for(int k=0;k<temp_M;k++){
                    auto subres=subpro(i,j,ii,k+1);
                    res[0][i][j][index_g++]=subres[0];
                    res[1][i][j][index_f++]=subres[1];
                    res[2][i][j][index_p++]=subres[2];
                }
            }
        }
    }

    return res;

}

vector<double> subpro(int area,int i,int area1,int j){
    vector<double> res(3,0);//G,f,p
    res[0]=D_MAX;
    if(j==-1){
        double f_local_min=max(f_range_local[area][i].first,D[area][i][1]*D[area][i][3]/D[area][i][4]);
        double f_local_max=min(f_range_local[area][i].second,sqrt(D[area][i][5]/(kn[area][i]*D[area][i][1]*D[area][i][3])));
        if(f_local_min<=f_local_max){
            double f_local_opt=pow(w_t[area][i]/(2*w_e[area][i]*kn[area][i]),1.0/3.0);
            
            if(f_local_opt<f_local_min){
                res[0]=func_local(area,i,f_local_min);
                res[1]=f_local_min;
            }
            else if(f_local_opt>f_local_max){
                res[0]=func_local(area,i,f_local_max);
                res[1]=f_local_max;
            }
            else{
                res[0]=func_local(area,i,f_local_opt);
                res[1]=f_local_opt;
            }
        }
    }
    else{
        if(area==area1){
            if(j==0){
                auto temp_res=f_to_x_opt_MBS(area,i);
                
                if(temp_res[0]!=0){
                    res[0]= func_MBS_SBS(area,i,j,temp_res[1],f_range_mec[area][j].first);
                   
                    res[2]= func_R_rev(area,i,j,temp_res[1]);
                    
                }
            }
            else{
                double x_SBS_min = func_R(area,i, j, p_range_all[area][i].second);
                double x_SBS_max = func_R(area,i, j, p_range_all[area][i].first);
                auto temp_res = f_to_x_range(area,i, j, x_SBS_min, x_SBS_max);
                if(temp_res[0]!=0){
                    auto temp_res1=fibonacci_search(area,i, j,temp_res[1], temp_res[2],fibonacci_itr);
                    res[0]=temp_res1[2];
                    res[1]=temp_res1[0];
                    res[2]=func_R_rev(area,i, j,temp_res1[1]);
                }
            }
        }
        else{
            if(j==0){

                auto temp_res = f_to_x_opt_MBS_1(area,i,area1);
                if(temp_res[0]!=0){
                    res[0] = func_MBS_SBS_1(area,i,area1,j,temp_res[1], f_range_mec[area1][j].first);
                    res[2] = temp_res[1];
                }
            }
            else{
                double p_SBS_min = p_range_all[area][i].first;
                double p_SBS_max = p_range_all[area][i].second;
                auto temp_res = f_to_x_range_1(area,i, area1,j, p_SBS_min, p_SBS_max);
                if(temp_res[0]!=0){
                    auto temp_res1=fibonacci_search_1(area,i,area1, j,temp_res[1], temp_res[2], fibonacci_itr);
                    // if(temp_res1[2]>D_MAX){
                    //     cout<<area<<i<<area1<<j<<endl;
                    // }
                    res[0]=temp_res1[2];
                    res[1]=temp_res1[0];
                    res[2]=temp_res1[1];
                }
            }
           
            
        }
    }
    return res;
    

}

//online
vector<vector<double>> exploration(vector<double> Nmax,int K,vector<int> KU,int SM,int Nm, double &Res,double &Res2){
    vector<vector<pair<int,int>>> memo_u(K);
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    auto temp_1=uniform2_int(0,SM-1,Area,Nm);
   for(int i=0;i<Area;i++){
        if(N[i]==0){
            continue;
        }
        res[i].resize(N[i]);
        res1[i].resize(N[i]);
        for(int j=0;j<Nmax[i];j++){
            int temp_j=temp_1[i][j];
            double temp_s=0;
            int index_j=0;
            for(int jj=0;jj<K;jj++){
                temp_s+=KU[jj];
                if(temp_s>temp_j){
                    index_j=jj;
                    memo_u[index_j].push_back({i,j});
                    break;
                } 
            }
                
        }
    }

    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        if(l1<=KU[j]){
            continue;
        }
        while(l1>KU[j]){
            //随机选择
            srand(seed);
            int index_temp=rand()%(l1);
            seed+=1;
            int index_a=memo_u[j][index_temp].first;
            int index_u=memo_u[j][index_temp].second;

            srand(seed);
            int index_temp_1=rand()%(Area+1);
            seed+=1;
           
  
            if(index_temp_1==0){
                memo_u[0].push_back({index_a,index_u});
            }
            else{
                index_temp_1=accumulate(ML.begin(),ML.begin()+index_temp_1-1,1);
                
                memo_u[index_temp_1].push_back({index_a,index_u});
            }
            memo_u[j].erase(memo_u[j].begin()+index_temp);
            l1--;
        }
    }
   // Print(ML);
    vector<double> temp_res;
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        for(int jj=0;jj<l1;jj++){
            int temp_a=memo_u[j][jj].first;
            int temp_u=memo_u[j][jj].second;
            int temp=N[temp_a];
            int temp_uu=temp_u%temp;
            
            if(j==0){
                
                temp_res=subpro(temp_a,temp_uu,temp_a,-1);
            }
            else{
            int temp_m=0;
            for(int area1=0;area1<Area;area1++){
                temp_m+=ML[area1];
                if(temp_m>j-1){
                    //cout<<j<<' '<<area1<<' '<<j-1-temp_m+ML[area1]<<endl;
                    temp_res=subpro(temp_a,temp_uu,area1,j-1-temp_m+ML[area1]);
                    break;
                }
            }  
            }

            V[temp_a][temp_u][j]+=1;
            if(v-p*log10(1+temp_res[0])<=0){
                cout<<"error"<<endl;
            }
            S[temp_a][temp_u][j]+=v-p*log10(1+temp_res[0]);
            if(temp_u<N[temp_a]){
                res1[temp_a][temp_u]=temp_res[0];
                res[temp_a][temp_u]=j;
            }
        }
    }
    //cout<<"e1"<<endl;
    for(auto i:res1){
        for(auto j:i){
           // cout<<j<<' ';
            Res2+=v-p*log10(1+j);
            Res+=j;
        }
    //cout<<endl;
     }
    // cout<<Res<<endl;
    return res;
}
vector<vector<double>> matching(double e,vector<int> KU,int SM,double& Res,double &Res2){
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    for(int area=0;area<Area;area++){
        if(N[area]==0){
            continue;
        }
        res[area].resize(N[area],-1);
        res1[area].resize(N[area]);
        for(int i=0;i<N[area];i++){
            if(a[area][i]==-1){
                double mymax=0,submax=0;
                vector<int> index_max;
                for(int m=0;m<SM;m++){
                    if(R1[area][i][m]-B[area][i][m]>mymax){
                        mymax=R1[area][i][m]-B[area][i][m];
                    }
                }
                for(int m=0;m<SM;m++){
                    if(R1[area][i][m]-B[area][i][m]==mymax){
                        index_max.push_back(m);
                    }
                    else{
                        if(R1[area][i][m]-B[area][i][m]>submax){
                            submax=R1[area][i][m]-B[area][i][m];
                        }
                    }
                }
                for(auto m:index_max){
                    B[area][i][m]=R1[area][i][m]-submax+e;
                }
                
            }
        }
    }
    int index_j=0;
    double temp_s=KU[0];
    vector<double> temp_res;
    for(int m=0;m<SM;m++){
        double val_max=0;
        int index_max_i=0;
        int index_area=0;
        for(int area=0;area<Area;area++){
            for(int i=0;i<N[area];i++){
                if(B[area][i][m]>val_max){
                    val_max=B[area][i][m];
                    index_max_i=i;
                    index_area=area;

                }
            }

        }
        if(temp_s==m){
            index_j=index_j+1;
            temp_s=temp_s+KU[index_j];
        }
        if(val_max!=0){
            if(index_j==0){
                temp_res=subpro(index_area,index_max_i,index_area,-1);
                a[index_area][index_max_i]=m;
                for(int area=0;area<Area;area++){
                    for(int i=0;i<N[area];i++){
                        if(a[area][i]==m&&area!=index_area&&i!=index_max_i){

                           a[area][i]=-1;
                        }       
                    }   
                }
                res1[index_area][index_max_i]=temp_res[0];
                res[index_area][index_max_i]=index_j;
            }
            else{
                double temp_m=0;
                for(int area1=0;area1<Area;area1++){
                    temp_m=temp_m+ML[area1];
                    if(temp_m>index_j-1){
                        temp_res=subpro(index_area,index_max_i,area1,index_j-1-temp_m+ML[area1]);
                        a[index_area][index_max_i]=m;
                        for(int area=0;area<Area;area++){
                            for(int i=0;i<N[area];i++){
                                if(a[area][i]==m&&area!=index_area&&i!=index_max_i){
                                    a[area][i]=-1;
                                }       
                            }   
                        }
                        res1[index_area][index_max_i]=temp_res[0];
                        
                        res[index_area][index_max_i]=index_j;
                        break;
                    }


                }
            }
        }
    }
    for(int area=0;area<Area;area++){
        for(int i=0;i<N[area];i++){
            if(res1[area][i]==0){
                //该用户被随机卸载
                srand(seed);
                int index_temp_1=rand()%(Area+1);
                seed+=1;
                if(index_temp_1==0){
                    temp_res=subpro(area,i,area,-1);
                    res1[area][i]=temp_res[0];
                    res[area][i]=0;
                }
                else{
                    temp_res=subpro(area,i,index_temp_1-1,0);
                    res1[area][i]=temp_res[0];
                    res[area][i]=accumulate(ML.begin(),ML.begin()+index_temp_1-1,1);
                
                }

            }
        }
    }
    // for(auto i:a){
    //     for(auto j:i){
    //         cout<<j<<' ';
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // Res=0;
   // cout<<"m"<<endl;
    for(auto i:res1){
        for(auto j:i){
           // cout<<j<<' ';
            Res2+=v-p*log10(1+j);
            Res+=j;
        }
   // cout<<endl;
     }
    // cout<<Res<<endl;
    return res;
}

vector<vector<double>> exploitation(int K,vector<int> temp_M,double& Res,double &Res2){
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    for(int area=0;area<Area;area++){
        if(N[area]==0){
            continue;
        }
        res[area].resize(N[area]);
        res1[area].resize(N[area]);
        for(int i=0;i<N[area];i++){
            auto temp_res=subpro(area,i,a1[area][i][0],a1[area][i][1]);
            res1[area][i]=temp_res[0];
            int temp=a[area][i];
            int temp_s=0;
            int index_j=0;
            for(int k=0;k<K;k++){
                temp_s+=temp_M[k];
                if(temp_s>temp){
                    index_j=k;
                    res[area][i]=index_j;
                    break;
                }
            }
            
        }
    } 
   // cout<<"e2"<<endl;
    for(auto i:res1){
        for(auto j:i){
           // cout<<j<<' ';
            Res+=j;
            Res2+=v-p*log10(1+j);
        }
       // cout<<endl;
    }
    //cout<<Res<<endl;
    return res;
}


//MUCB
vector<vector<double>>  MUCB(double t,vector<double> Nmax,int K,vector<int> temp_M,double &Res,double &Res2){
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    vector<vector<pair<int,int>>> memo_u(K);
    //先用户选择,每次做最大数量用户
    for(int area=0;area<Area;area++){
        if(N[area]==0){
            continue;
        }
        res[area].resize(N[area]);
        res1[area].resize(N[area]);
        for(int i=0;i<Nmax[area];i++){
            double temp_max=0.0;
            int temp_index=0;
            for(int j=0;j<K;j++){
                if(C_V[area][i][j][0]==0){
                    temp_index=j;
                    break;
                }
                else{
                    double temp_p=C_V[area][i][j][1]+sqrt(2.0*log(t)/C_V[area][i][j][0]);
                    if(temp_max<temp_p){
                        temp_max=temp_p;
                        temp_index=j;
                    }
                }
            }
            memo_u[temp_index].push_back({area,i});
        }
    }
    //超过上限的用户再随机选择
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        if(l1<=temp_M[j]){
            continue;
        }
        while(l1>temp_M[j]){
            //随机选择
            srand(seed);
            int index_temp=rand()%(l1);
            seed+=1;
            int index_a=memo_u[j][index_temp].first;
            int index_u=memo_u[j][index_temp].second;

            srand(seed);
            int index_temp_1=rand()%(Area+1);
            seed+=1;
            if(index_temp_1==0){
                memo_u[0].push_back({index_a,index_u});
            }
            else{
                index_temp_1=accumulate(ML.begin(),ML.begin()+index_temp_1-1,1);
                
                memo_u[index_temp_1].push_back({index_a,index_u});
            }
            memo_u[j].erase(memo_u[j].begin()+index_temp);
            l1--;
        }
    }
    //服务器决定
    vector<double> temp_res;
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        for(int jj=0;jj<l1;jj++){
            int temp_a=memo_u[j][jj].first;
            int temp_u=memo_u[j][jj].second;
            int temp=N[temp_a];
            int temp_uu=temp_u%temp;  
            if(j==0){
                temp_res=subpro(temp_a,temp_uu,temp_a,-1);
            }
            else{
            int temp_m=0;
            for(int area1=0;area1<Area;area1++){
                temp_m+=ML[area1];
                if(temp_m>j-1){
                    //cout<<j<<' '<<area1<<' '<<j-1-temp_m+ML[area1]<<endl;
                    temp_res=subpro(temp_a,temp_uu,area1,j-1-temp_m+ML[area1]);
                    break;
                }
            }  
            }

            C_V[temp_a][temp_u][j][0]+=1;
            if(v-p*log10(1+temp_res[0])<=0){
                cout<<"error"<<endl;
            }
            C_V[temp_a][temp_u][j][1]=(t-1)*C_V[temp_a][temp_u][j][1]/t+(v-p*log10(1+temp_res[0]))/(v*t);
            if(temp_u<N[temp_a]){
                res1[temp_a][temp_u]=temp_res[0];
                res[temp_a][temp_u]=j;
            }
        }
    }
    for(auto i:res1){
        for(auto j:i){
            Res+=j;
            Res2+=v-p*log10(1+j);
        }
       
    }
    return res;



}
 
//MEXP3
vector<vector<double>>  MEXP3(double t,vector<double> Nmax,int Nm,int K,vector<int> temp_M,double &Res,double &Res2){
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    vector<vector<pair<int,int>>> memo_u(K);
    vector<vector<vector<double>>> P(Area);
    double h=sqrt(log(K)/(t*K));
    //先用户选择,每次做最大数量用户
    //随机数
    auto temp_rand=uniform2(0.0,1.0,Area,Nm);
    for(int area=0;area<Area;area++){
        if(N[area]==0){
            continue;
        }
        res[area].resize(N[area]);
        res1[area].resize(N[area]);
        P[area].resize(Nmax[area],vector<double>(K,0));
        for(int i=0;i<Nmax[area];i++){
            //计算该时刻该用户的P
            double temp=0;
            for(int j=0;j<K;j++){
                
                temp+=exp(h*S1[area][i][j]);
            }
            double temp_1=0.0;
            for(int j=0;j<K;j++){
                P[area][i][j]=exp(h*S1[area][i][j])/temp;
                if(temp_rand[area][i]>=temp_1&&temp_rand[area][i]<temp_1+P[area][i][j]){
                memo_u[j].push_back({area,i});
                break;
               }
               temp_1+=P[area][i][j];
               
            }
            
        }
    }
    //超过上限的用户再随机选择
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        if(l1<=temp_M[j]){
            continue;
        }
        while(l1>temp_M[j]){
            //随机选择
            srand(seed);
            int index_temp=rand()%(l1);
            seed+=1;
            int index_a=memo_u[j][index_temp].first;
            int index_u=memo_u[j][index_temp].second;

            srand(seed);
            int index_temp_1=rand()%(Area+1);
            seed+=1;
            if(index_temp_1==0){
                memo_u[0].push_back({index_a,index_u});
            }
            else{
                index_temp_1=accumulate(ML.begin(),ML.begin()+index_temp_1-1,1);
                
                memo_u[index_temp_1].push_back({index_a,index_u});
            }
            memo_u[j].erase(memo_u[j].begin()+index_temp);
            l1--;
        }
    }
    for(int area=0;area<Area;area++){
        for(int i=0;i<Nmax[area];i++){
            for(int j=0;j<K;j++){
                S1[area][i][j]+=1;
            }
        }
    }
    //服务器决定
    vector<double> temp_res;
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        for(int jj=0;jj<l1;jj++){
            int temp_a=memo_u[j][jj].first;
            int temp_u=memo_u[j][jj].second;
            int temp=N[temp_a];
            int temp_uu=temp_u%temp;  
            if(j==0){
                temp_res=subpro(temp_a,temp_uu,temp_a,-1);
            }
            else{
            int temp_m=0;
            for(int area1=0;area1<Area;area1++){
                temp_m+=ML[area1];
                if(temp_m>j-1){
                    //cout<<j<<' '<<area1<<' '<<j-1-temp_m+ML[area1]<<endl;
                    temp_res=subpro(temp_a,temp_uu,area1,j-1-temp_m+ML[area1]);
                    break;
                }
            }  
            }

            
            if((v-p*log10(1+temp_res[0]))/v<=0||(v-p*log10(1+temp_res[0]))/v>1){
                cout<<"error"<<endl;
            }
            S1[temp_a][temp_u][j]=S1[temp_a][temp_u][j]-(1-(v-p*log10(1+temp_res[0]))/v)/P[temp_a][temp_u][j];
            if(temp_u<N[temp_a]){
                res1[temp_a][temp_u]=temp_res[0];
                res[temp_a][temp_u]=j;
            }
        }
    }
    for(auto i:res1){
        for(auto j:i){
            Res+=j;
            Res2+=v-p*log10(1+j);
        }
       
    }
    return res;



}


//MEXP3-iX
vector<vector<double>>  MEXP3_ix(double t,vector<double> Nmax,int Nm,int K,vector<int> temp_M,double &Res,double &Res2){
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    vector<vector<pair<int,int>>> memo_u(K);
    vector<vector<vector<double>>> P(Area);
    double h=sqrt(2.0*log(K+1)/(t*K));
    double r=h/2.0;
    //先用户选择,每次做最大数量用户
    //随机数
    auto temp_rand=uniform2(0.0,1.0,Area,Nm);
    for(int area=0;area<Area;area++){
        if(N[area]==0){
            continue;
        }
        res[area].resize(N[area]);
        res1[area].resize(N[area]);
        P[area].resize(Nmax[area],vector<double>(K,0));
        for(int i=0;i<Nmax[area];i++){
            //计算该时刻该用户的P
            double temp=0;
            for(int j=0;j<K;j++){                
                temp+=exp(-h*S2[area][i][j]);
            }
            double temp_1=0.0;
            for(int j=0;j<K;j++){
                P[area][i][j]=exp(-h*S2[area][i][j])/temp;
                if(temp_rand[area][i]>=temp_1&&temp_rand[area][i]<temp_1+P[area][i][j]){
                memo_u[j].push_back({area,i});
                break;
               }
               temp_1+=P[area][i][j];
               
            }
            
        }
    }
    //超过上限的用户再随机选择
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        if(l1<=temp_M[j]){
            continue;
        }
        while(l1>temp_M[j]){
            //随机选择
            srand(seed);
            int index_temp=rand()%(l1);
            seed+=1;
            int index_a=memo_u[j][index_temp].first;
            int index_u=memo_u[j][index_temp].second;

            srand(seed);
            int index_temp_1=rand()%(Area+1);
            seed+=1;
            if(index_temp_1==0){
                memo_u[0].push_back({index_a,index_u});
            }
            else{
                index_temp_1=accumulate(ML.begin(),ML.begin()+index_temp_1-1,1);
                
                memo_u[index_temp_1].push_back({index_a,index_u});
            }
            memo_u[j].erase(memo_u[j].begin()+index_temp);
            l1--;
        }
    }
    //服务器决定
    vector<double> temp_res;
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        for(int jj=0;jj<l1;jj++){
            int temp_a=memo_u[j][jj].first;
            int temp_u=memo_u[j][jj].second;
            int temp=N[temp_a];
            int temp_uu=temp_u%temp;  
            if(j==0){
                temp_res=subpro(temp_a,temp_uu,temp_a,-1);
            }
            else{
            int temp_m=0;
            for(int area1=0;area1<Area;area1++){
                temp_m+=ML[area1];
                if(temp_m>j-1){
                    //cout<<j<<' '<<area1<<' '<<j-1-temp_m+ML[area1]<<endl;
                    temp_res=subpro(temp_a,temp_uu,area1,j-1-temp_m+ML[area1]);
                    break;
                }
            }  
            }

            S2[temp_a][temp_u][j]+=(1-(v-p*log10(1+temp_res[0]))/v)/(P[temp_a][temp_u][j]+r);
            if(temp_u<N[temp_a]){
                res1[temp_a][temp_u]=temp_res[0];
                res[temp_a][temp_u]=j;
            }
        }
    }
    for(auto i:res1){
        for(auto j:i){
            Res+=j;
            Res2+=v-p*log10(1+j);
        }
       
    }
    return res;



}



double JUACO(vector<int> temp_M){
    int M_s=accumulate(M.begin(),M.end(),0.0);
    int N_s=accumulate(N.begin(),N.end(),0.0);
    double res1=0.0,res2=1.0;
    double e1=0.00001,e2=0.00001;
    vector<vector<vector<double>>> x(Area);
    vector<vector<vector<double>>> f(Area);
    vector<vector<vector<double>>> p(Area);
    for(int i=0;i<Area;i++){
        if(N[i]!=0){
            x[i].resize(N[i],vector<double>(M_s+Area+1));
            f[i].resize(N[i],vector<double>(M_s+1));
            p[i].resize(N[i],vector<double>(M_s+Area));
            for(int j=0;j<N[i];j++){
                f[i][j][0]=f_range_local[i][j].first;
                for(int ii=0,index=1;ii<Area;ii++){
                    for(int jj=0;jj<M[ii];jj++){
                        f[i][j][index++]=f_range_mec[ii][jj+1].first;
                    }
                }
                for(int ii=0;ii<Area+M_s;ii++){
                    p[i][j][ii]=p_range_all[i][j].first;
                }

            }
        }
    }
    
    vector<vector<double>> G(N_s,vector<double>(M_s+Area+1));
    do{
        //x
        res1=res2;
        for(int i=0,index=0;i<Area;i++){
            for(int j=0;j<N[i];j++,index++){
                G[index][0]=func_local(i,j,f[i][j][0]);
                int index1=1;
                int index2=1;
                for(int ii=0;ii<Area;ii++){
                    G[index][index1]=func_MBS_SBS_11(i,j,ii,0,p[i][j][index1-1],f_range_mec[ii][0].first);
                    index1++;
                    for(int jj=1;jj<=M[ii];jj++){
                        G[index][index1]=func_MBS_SBS_11(i,j,ii,jj,p[i][j][index1-1],f[i][j][index2++]);
                        index1++;
                    }
                }
            }
        }
        double res3=0.0,res4=0.0;
        vector<vector<double>> G_1(N_s);
        for(int j=0;j<N_s;j++){
            for(int jj=0;jj<G[0].size();jj++){
                for(int k=0;k<temp_M[jj];k++){
                    G_1[j].push_back(G[j][jj]);
                }
            }
        }
        lenm=G_1[0].size();
        lenn=G_1.size();
        for (int j = 0; j < lenn; ++j){
            for (int jj = 0; jj < lenm; ++jj){
                E_val[j][jj]=10000-G_1[j][jj]; 
                if(E_val[j][jj]<0){
                    cout<<"error"<<endl;
                }    
            }
                    
        }
        for (int j = lenn; j < lenm; ++j){
            for (int jj = 0; jj < lenm; ++jj){
                E_val[j][jj]=0;
            }
                    
        }
        double Opt=0;
        auto opt_s=KM(Opt);
        res4=10000*lenn-Opt;
        vector<vector<vector<int>>> U(Area);
        for(int i=0,index=0;i<Area;i++){
            U[i].resize(N[i],vector<int>(2,-2));
            for(int j=0;j<N[i];j++){
                int temp_j=0;
                for(int k=0;k<lenm;k++){
                    if(opt_s[index][k]==1){
                        temp_j=k;
                        break;
                    }
                }
                index++;
                double temp_s=0;
                int index_j=0;
                for(int jj=0;jj<temp_M.size();jj++){
                    temp_s+=temp_M[jj];
                    if(temp_s>temp_j){
                        index_j=jj;
                        break;
                    } 
                }
                if(index_j==0){
                    //local
                    U[i][j][0]=i;
                    U[i][j][1]=-1;
                        
                }
                else{
                    int temp_m=0;
                    for(int area1=0;area1<Area;area1++){
                        temp_m+=ML[area1];
                        if(temp_m>index_j-1){ 
                            //area1 cell index_j-1-temp_m+ML[area1] 服务器 
                            U[i][j][0]=area1;
                            U[i][j][1]=index_j-1-temp_m+ML[area1];  
                            break;
                        }
                    } 

                }
            }
        }
        //f,p
        do{
            res3=res4;
            res4=0;//重新优化所有的用户决策目标
            for(int area=0;area<Area;area++){
                for(int i=0;i<N[area];i++){
                    int area1=U[area][i][0];
                    int j=U[area][i][1];
                    //f
                    if(j==-1){
                        //本地
                        double f_local_min=max(f_range_local[area][i].first,D[area][i][1]*D[area][i][3]/D[area][i][4]);
                        double f_local_max=min(f_range_local[area][i].second,sqrt(D[area][i][5]/(kn[area][i]*D[area][i][1]*D[area][i][3])));
                        if(f_local_min<=f_local_max){
                            double f_local_opt=pow(w_t[area][i]/(2*w_e[area][i]*kn[area][i]),1.0/3.0);
                            if(f_local_opt<f_local_min){
                                res4+=func_local(area,i,f_local_min);
                                f[area][i][0]=f_local_min;
                            }
                            else if(f_local_opt>f_local_max){
                                res4+=func_local(area,i,f_local_max);
                                f[area][i][0]=f_local_max;
                            }
                            else{
                                res4+=func_local(area,i,f_local_opt);
                                f[area][i][0]=f_local_opt;
                            }
                        }

                    }
                    else if(j==0){
                        //不需要优化f
                    }
                    else{
                        //边缘服务器
                        int temp_j=accumulate(M.begin(),M.begin()+area1,j);//对于服务器的f下标
                        int temp_jj=accumulate(ML.begin(),ML.end()+area1,j);//对于服务器的p下标
                        double p1=p[area][i][temp_jj];
                        if(area==area1){
                            double f_local_min=max(f_range_mec[area][j].first,D[area][i][1]*D[area][i][3]/(D[area][i][4]-(D[area][i][0]+D[area][i][2])*func_R(area,i,j,p1)));
                            double f_local_max=min(f_range_mec[area][j].second,sqrt((D[area][i][5]-(D[area][i][0]+D[area][i][2])*func_R(area,i,j,p1)*p1)/(km[area][j]*D[area][i][1]*D[area][i][3])));
                            if(f_local_min<=f_local_max){
                                double f_local_opt=pow((w_t[area][i]+w_c[area][i]*cost[area][j])/(2*w_e[area][i]*km[area][j]),1.0/3.0);
                                if(f_local_opt<f_local_min){
                                    
                                    f[area][i][temp_j]=f_local_min;
                                }  
                                else if(f_local_opt>f_local_max){
                                    f[area][i][temp_j]=f_local_max;
                                }
                                else{
                                    f[area][i][temp_j]=f_local_opt;
                                }
                            }
                            else{
                                cout<<"e1"<<endl;
                            }
                        }
                        else{//跨单元
                            double f_local_min=max(f_range_mec[area1][j].first,D[area][i][1]*D[area][i][3]/(D[area][i][4]-(D[area][i][0]+D[area][i][2])*func_R_1(area,i,area1,j,p1)));
                            double f_local_max=min(f_range_mec[area1][j].second,sqrt((D[area][i][5]-(D[area][i][0]+D[area][i][2])*func_R_1(area,i,area1,j,p1)*p1)/(km[area1][j]*D[area][i][1]*D[area][i][3])));
                            if(f_local_min<=f_local_max){
                                double f_local_opt=pow((w_t[area][i]+w_c[area][i]*cost[area1][j])/(2*w_e[area][i]*km[area1][j]),1.0/3.0);
                                if(f_local_opt<f_local_min){
                                    
                                    f[area][i][temp_j]=f_local_min;
                                }  
                                else if(f_local_opt>f_local_max){
                                    f[area][i][temp_j]=f_local_max;
                                }
                                else{
                                    f[area][i][temp_j]=f_local_opt;
                                }
                            }
                            else{
                                cout<<"e2"<<endl;
                            }

                        }
                    }
                    //p
                     
                    if(j!=-1){
                        if(j==0){
                            //云服务器
                            int temp_j=accumulate(ML.begin(),ML.begin()+area1,0);
                            if(area==area1){
                                auto temp_res=f_to_x_opt_MBS(area,i);
                                if(temp_res[0]!=0){
                                    p[area][i][temp_j]= func_R_rev(area,i,j,temp_res[1]);
                                    res4+=func_MBS_SBS_11(area,i,area1,j,p[area][i][temp_j],f_range_mec[area1][j].first);
                                }
                                else{
                                    cout<<"e3"<<endl;
                                }
                            }
                            else{
                                auto temp_res = f_to_x_opt_MBS_1(area,i,area1);
                                if(temp_res[0]!=0){
                                    p[area][i][temp_j] = temp_res[1];
                                    res4+=func_MBS_SBS_11(area,i,area1,j,temp_res[1], f_range_mec[area1][j].first);
                                }
                                else{
                                    cout<<"e4"<<endl;
                                }
 
                            }
                            
                        }
                        else{
                            //边缘服务器,存在干扰
                            int temp_jj=accumulate(M.begin(),M.begin()+area1,j);//对于服务器的f下标
                            int temp_j=accumulate(ML.begin(),ML.begin()+area1,j);
                            double f1=f[area][i][temp_jj];
                            if(area==area1){
                                auto temp_res=f_to_x_opt_MBS1(area,i,j,f1);
                                if(temp_res[0]!=0){
                                    p[area][i][temp_j]= func_R_rev(area,i,j,temp_res[1]);
                                    res4+=func_MBS_SBS_11(area,i,area1,j,p[area][i][temp_j],f1);
                                }
                                else{
                                    cout<<"e5"<<endl;
                                }
                            }
                            else{
                                //无干扰
                                auto temp_res = f_to_x_opt_MBS1_1(area,i,area1,j,f1);
                                if(temp_res[0]!=0){
                                    p[area][i][temp_j] = temp_res[1];
                                    res4+=func_MBS_SBS_11(area,i,area1,j,temp_res[1], f1);
                                }
                                else{
                                    cout<<"e6"<<endl;
                                }
 
                            }
                        }
                    }
                    
                }
            }
        }while(abs(res4-res3)>e1);
        res2=res4;
    }while(abs((res2-res1)/res1)>e2);
    return res2;
}

double func_obj(int area,int i,int area1,int j,vector<double> x){//f,p
    if(j==-1){
        return func_local(area,i,x[0]);
    }
    else{
        if(j==0){
            return func_MBS_SBS_11(area,i,area1,j,x[0],f_range_mec[area1][j].first);
        }
        else{
            return func_MBS_SBS_11(area,i,area1,j,x[1],x[0]);
        }
    }

}

bool  conl_T(int area,int i,int area1,int j,vector<double> x){
    if(j==-1){
        if(D[area][i][1]*D[area][i][3]/x[0]<=D[area][i][4]){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        if(area==area1){
            if(j==0){
                if((D[area][i][0]+D[area][i][2])*func_R(area,i,j,x[0])+D[area][i][1]*D[area][i][3]/f_range_mec[area][j].first<=D[area][i][4]){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                if((D[area][i][0]+D[area][i][2])*func_R(area,i,j,x[1])+D[area][i][1]*D[area][i][3]/x[0]<=D[area][i][4]){
                    return true;
                }
                else{
                    return false;
                }
            }

        }
        else{
           if(j==0){
                if((D[area][i][0]+D[area][i][2])*func_R_1(area,i,area1,j,x[0])+D[area][i][1]*D[area][i][3]/f_range_mec[area1][j].first<=D[area][i][4]){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                if((D[area][i][0]+D[area][i][2])*func_R_1(area,i,area1,j,x[1])+D[area][i][1]*D[area][i][3]/x[0]<=D[area][i][4]){
                    return true;
                }
                else{
                    return false;
                }
            }

        }


    }


}
bool  conl_E(int area,int i,int area1,int j,vector<double> x){
    if(j==-1){
        if(kn[area][i]*D[area][i][1]*D[area][i][3]*x[0]*x[0]<=D[area][i][5]){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        if(area==area1){
            if(j==0){
                if((D[area][i][0]+D[area][i][2])*func_R(area,i,j,x[0])*x[0]+km[area][j]*D[area][i][1]*D[area][i][3]*f_range_mec[area1][j].first*f_range_mec[area1][j].first<=D[area][i][5]){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                if((D[area][i][0]+D[area][i][2])*func_R(area,i,j,x[1])*x[1]+km[area][j]*D[area][i][1]*D[area][i][3]*x[0]*x[0]<=D[area][i][5]){
                    return true;
                }
                else{
                    return false;
                }
            }

        }
        else{
           if(j==0){
                if((D[area][i][0]+D[area][i][2])*func_R_1(area,i,area1,j,x[0])*x[0]+km[area1][j]*D[area][i][1]*D[area][i][3]*f_range_mec[area1][j].first*f_range_mec[area1][j].first<=D[area][i][5]){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                if((D[area][i][0]+D[area][i][2])*func_R_1(area,i,area1,j,x[1])*x[1]+km[area1][j]*D[area][i][1]*D[area][i][3]*x[0]*x[0]<=D[area][i][5]){
                    return true;
                }
                else{
                    return false;
                }
            }

        }


    }

}



vector<double> pso(int area,int i,int area1,int j,int dim,vector<double> xlimit_min,vector<double> xlimit_max,bool (*func1)(int,int,int,int,vector<double>),bool (*func2)(int,int,int,int,vector<double>),double (*func)(int,int,int,int,vector<double>)){
    vector<double> res(dim+1);
    res[0]=D_MAX;
    vector<double> vlimit_max(dim,1);
    vector<double> vlimit_min(dim,-1);
    vector<vector<double>> pop_x(popsize,vector<double>(dim,0));
    vector<vector<double>> pop_v(popsize,vector<double>(dim,0));
    vector<double> fitness_pop(popsize,0);
    vector<double> fitness_lbest(popsize,0);
    for(int k=0;k<popsize;k++){
        // 位置初始化
        for(int kk=0;kk<dim;kk++){
            pop_x[k][kk] = xlimit_min[kk] + prr[0][k][kk]*(xlimit_max[kk] - xlimit_min[kk]);
        }
        // 速度初始化
        for(int kk=0;kk<dim;kk++){
            pop_v[k][kk] = vlimit_min[kk] + prr[0][k][kk]*(vlimit_max[kk] - vlimit_min[kk]);
        }
    }

    // 初始化个体极值
    vector<vector<double>> lbest(pop_x); // 个体历史最佳极值记录
    for(int k=0;k<popsize;k++){ 
        if((*func1)(area,i,area1,j,pop_x[k])){
            if((*func2)(area,i,area1,j,pop_x[k])){
                fitness_lbest[k] = (*func)(area,i,area1,j,pop_x[k]);
            }
            else{
                fitness_lbest[k] = D_MAX;
            }
        }
        else{
            fitness_lbest[k] = D_MAX;
        }      
    }

    // 初始化全局极值
    vector<double> popbest(pop_x[0]);
    double fitness_popbest = fitness_lbest[0];
    for(int k=1;k<popsize;k++){ 
        if(fitness_lbest[k] < fitness_popbest){
            fitness_popbest = fitness_lbest[k];
            popbest = pop_x[k];
        }
    }

    // 粒子群迭代
    vector<double> record; // 记录每次迭代的全局极小值
    for(int iter=0;iter<max_iter;iter++){
        for(int k=0;k<popsize;k++){ 
            // 更新速度 边界处理
            for(int kk=0;kk<dim;kk++){
                pop_v[k][kk] = w1*pop_v[k][kk] + c1*prr[iter][k][kk]*(lbest[k][kk] - pop_x[k][kk]) + c2*prr[iter][k][kk]*(popbest[kk] - pop_x[k][kk]);
            }
            
            for(int kk=0;kk<dim;kk++){ 
                if(pop_v[k][kk] > vlimit_max[kk]) 
                    pop_v[k][kk] = vlimit_max[kk];
                else if(pop_v[k][kk] < vlimit_min[kk]){
                    pop_v[k][kk] = vlimit_min[kk];
                }
            }
            // 更新位置 边界处理 修正位置 （等式约束）
            for(int kk=0;kk<dim;kk++){
                pop_x[k][kk] = pop_x[k][kk]+pop_v[k][kk];
            }
            
            for(int kk=0;kk<dim;kk++){ 
                if(pop_x[k][kk] > xlimit_max[kk]) 
                    pop_x[k][kk] = xlimit_max[kk];
                else if(pop_x[k][kk] < xlimit_min[kk]){
                    pop_x[k][kk] = xlimit_min[kk];
                }
            }  
            if((*func1)(area,i,area1,j,pop_x[k])){
                if((*func2)(area,i,area1,j,pop_x[k])){
                    fitness_pop[k] = (*func)(area,i,area1,j,pop_x[k]);
                }
                else{
                    fitness_pop[k] = D_MAX;
                }
            }
            else{
                fitness_pop[k] = D_MAX;
            }      
            // 当前适应度与个体历史最佳适应度作比较
            if(fitness_pop[k] < fitness_lbest[k]){
                for(int kk=0;kk<dim;kk++){
                    lbest[k][kk] = pop_x[k][kk];
                }
                fitness_lbest[k] = fitness_pop[k];
            }
            // 个体历史最佳适应度与种群历史最佳适应度作比较
            if(fitness_popbest > fitness_lbest[k]){
                fitness_popbest = fitness_lbest[k];
                for(int kk=0;kk<dim;kk++){
                    popbest[kk] = lbest[k][kk];
                }

            }
        }
        record.push_back(fitness_popbest);
    }
    res[0]= fitness_popbest;
    for(int k=0;k<dim;k++){
        res[k+1]=popbest[k];
    }
    return res;
}

void pso_thread(int area,int i){
    //每个用户执行一个线程
    int M_s=accumulate(M.begin(),M.end(),0);
    vector<double> res_g(M_s+Area+1);//G
    vector<double> res_f(M_s+1);//f
    vector<double> res_p(M_s+Area);//p
    int index_g=0,index_p=0,index_f=0;  
    auto subres=pso(area,i,area,-1,1,{f_range_local[area][i].first},{f_range_local[area][i].second},conl_T,conl_E,func_obj);
    res_g[index_g++]=subres[0];
    res_f[index_f++]=subres[1];
    for(int ii=0;ii<Area;ii++){
        int temp_M=M[ii];
        auto subres=pso(area,i,ii,0,1,{p_range_all[area][i].first},{p_range_all[area][i].second},conl_T,conl_E,func_obj);
        res_g[index_g++]=subres[0];
        res_p[index_p++]=subres[1];
        for(int k=0;k<temp_M;k++){
            if(F_MAX[ii][k]==0){
                index_g++;
                index_f++;
                index_p++;
            }
            else{
                auto subres=pso(area,i,ii,k+1,2,{f_range_mec[ii][k+1].first,p_range_all[area][i].first},{f_range_mec[ii][k+1].second,p_range_all[area][i].second},conl_T,conl_E,func_obj);
                res_g[index_g++]=subres[0];
                res_f[index_f++]=subres[1];
                res_p[index_p++]=subres[2];
            }
                    
        }
    }
    res111[0][area][i]=res_g;
    res111[1][area][i]=res_f;
    res111[2][area][i]=res_p;
}

vector<vector<vector<vector<double>>>> PSO(){
    //每个子问题用算法求解
    vector<vector<vector<vector<double>>>> res(3,vector<vector<vector<double>>>(Area));
    int M_s=accumulate(M.begin(),M.end(),0.0);
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        res[0][i].resize(temp_N,vector<double>(M_s+Area+1));//G
        res[1][i].resize(temp_N,vector<double>(M_s+1));//f
        res[2][i].resize(temp_N,vector<double>(M_s+Area));//p
        for(int j=0;j<temp_N;j++){
            int index_g=0,index_p=0,index_f=0;  
            auto subres=pso(i,j,i,-1,1,{f_range_local[i][j].first},{f_range_local[i][j].second},conl_T,conl_E,func_obj);
            res[0][i][j][index_g++]=subres[0];
            res[1][i][j][index_f++]=subres[1];
            for(int ii=0;ii<Area;ii++){
                int temp_M=M[ii];
                auto subres=pso(i,j,ii,0,1,{p_range_all[i][j].first},{p_range_all[i][j].second},conl_T,conl_E,func_obj);
                res[0][i][j][index_g++]=subres[0];
                res[2][i][j][index_p++]=subres[1];
                for(int k=0;k<temp_M;k++){
                    if(F_MAX[ii][k]==0){
                        index_g++;
                        index_f++;
                        index_p++;
                    }
                    else{
                        auto subres=pso(i,j,ii,k+1,2,{f_range_mec[ii][k+1].first,p_range_all[i][j].first},{f_range_mec[ii][k+1].second,p_range_all[i][j].second},conl_T,conl_E,func_obj);
                        res[0][i][j][index_g++]=subres[0];
                        res[1][i][j][index_f++]=subres[1];
                        res[2][i][j][index_p++]=subres[2];
                    }
                    
                }
            }
        }
    }

    return res;
}

void PSO_thread(){
    res111.clear();
    res111.resize(3,vector<vector<vector<double>>>(Area));
    int M_s=accumulate(M.begin(),M.end(),0.0);
    prr.resize(max_iter);
    for(int i=0;i<max_iter;i++){
        prr[i]=uniform2(0,1,popsize,2);
    }
    vector<thread> thread_vec; 
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        res111[0][i].resize(temp_N,vector<double>(M_s+Area+1));//G
        res111[1][i].resize(temp_N,vector<double>(M_s+1));//f
        res111[2][i].resize(temp_N,vector<double>(M_s+Area));//p
    }
    for(int i=0;i<Area;i++){
        int temp_N=N[i];
        if(temp_N==0){
            continue;
        }
        for(int j=0;j<temp_N;j++){
            thread_vec.emplace_back(thread(pso_thread,i,j));
        }
    }
    for (thread& t : thread_vec) 
       t.join(); 
}

//D-DEBO
vector<vector<double>> exploration1(int K,vector<int> KU,int SM, double &Res,double &Res2){
    vector<vector<pair<int,int>>> memo_u(K);
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    int Nm=0;
    for(int i=0;i<Area;i++){
        if(Nm<N[i]){
            Nm=N[i];
        }
    }
    auto temp_1=uniform2_int(0,SM-1,Area,Nm);
   for(int i=0;i<Area;i++){
        if(N[i]==0){
            continue;
        }
        res[i].resize(N[i]);
        res1[i].resize(N[i]);
        for(int j=0;j<N[i];j++){
            int temp_j=temp_1[i][j];
            double temp_s=0;
            int index_j=0;
            for(int jj=0;jj<K;jj++){
                temp_s+=KU[jj];
                if(temp_s>temp_j){
                    index_j=jj;
                    memo_u[index_j].push_back({i,j});
                    break;
                } 
            }
                
        }
    }

    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        if(l1<=KU[j]){
            continue;
        }
        while(l1>KU[j]){
            //随机选择
            srand(seed);
            int index_temp=rand()%(l1);
            seed+=1;
            int index_a=memo_u[j][index_temp].first;
            int index_u=memo_u[j][index_temp].second;

            srand(seed);
            int index_temp_1=rand()%(Area+1);
            seed+=1;
           
  
            if(index_temp_1==0){
                memo_u[0].push_back({index_a,index_u});
            }
            else{
                index_temp_1=accumulate(ML.begin(),ML.begin()+index_temp_1-1,1);
                
                memo_u[index_temp_1].push_back({index_a,index_u});
            }
            memo_u[j].erase(memo_u[j].begin()+index_temp);
            l1--;
        }
    }
   // Print(ML);
    vector<double> temp_res;
    for(int j=0;j<K;j++){
        int l1=memo_u[j].size();
        for(int jj=0;jj<l1;jj++){
            int temp_a=memo_u[j][jj].first;
            int temp_u=memo_u[j][jj].second;
            if(j==0){
                temp_res=subpro(temp_a,temp_u,temp_a,-1);
            }
            else{
            int temp_m=0;
            for(int area1=0;area1<Area;area1++){
                temp_m+=ML[area1];
                if(temp_m>j-1){
                    //cout<<j<<' '<<area1<<' '<<j-1-temp_m+ML[area1]<<endl;
                    temp_res=subpro(temp_a,temp_u,area1,j-1-temp_m+ML[area1]);
                    break;
                }
            }  
            }

            V_d[temp_a][temp_u][j]+=1;
            if(v-p*log10(1+temp_res[0])<=0){
                cout<<"error"<<endl;
            }
            S_d[temp_a][temp_u][j]+=v-p*log10(1+temp_res[0]);
            res1[temp_a][temp_u]=temp_res[0];
            res[temp_a][temp_u]=j;
        }
    }
    //cout<<"e11"<<endl;
    for(auto i:res1){
        for(auto j:i){
            //cout<<j<<' ';
            Res2+=v-p*log10(1+j);
            Res+=j;
        }
    //cout<<endl;
     }
     //cout<<Res<<endl;
    return res;
}
vector<vector<double>> matching1(double e,vector<int> KU,int SM,double& Res,double &Res2){
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    for(int area=0;area<Area;area++){
        if(N[area]==0){
            continue;
        }
        res[area].resize(N[area],-1);
        res1[area].resize(N[area]);
        for(int i=0;i<N[area];i++){
            if(a_d[area][i]==-1){
                double mymax=0,submax=0;
                vector<int> index_max;
                for(int m=0;m<SM;m++){
                    if(R1_d[area][i][m]-B_d[area][i][m]>mymax){
                        mymax=R1_d[area][i][m]-B_d[area][i][m];
                    }
                }
                for(int m=0;m<SM;m++){
                    if(R1_d[area][i][m]-B_d[area][i][m]==mymax){
                        index_max.push_back(m);
                    }
                    else{
                        if(R1_d[area][i][m]-B_d[area][i][m]>submax){
                            submax=R1_d[area][i][m]-B_d[area][i][m];
                        }
                    }
                }
                for(auto m:index_max){
                    B_d[area][i][m]=R1_d[area][i][m]-submax+e;
                }
                
            }
        }
    }
    int index_j=0;
    double temp_s=KU[0];
    vector<double> temp_res;
    for(int m=0;m<SM;m++){
        double val_max=0;
        int index_max_i=0;
        int index_area=0;
        for(int area=0;area<Area;area++){
            for(int i=0;i<N[area];i++){
                if(B_d[area][i][m]>val_max){
                    val_max=B_d[area][i][m];
                    index_max_i=i;
                    index_area=area;

                }
            }

        }
        if(temp_s==m){
            index_j=index_j+1;
            temp_s=temp_s+KU[index_j];
        }
        if(val_max!=0){
            if(index_j==0){
                temp_res=subpro(index_area,index_max_i,index_area,-1);
                a_d[index_area][index_max_i]=m;
                for(int area=0;area<Area;area++){
                    for(int i=0;i<N[area];i++){
                        if(a_d[area][i]==m&&area!=index_area&&i!=index_max_i){

                           a_d[area][i]=-1;
                        }       
                    }   
                }
                res1[index_area][index_max_i]=temp_res[0];
                res[index_area][index_max_i]=index_j;
            }
            else{
                double temp_m=0;
                for(int area1=0;area1<Area;area1++){
                    temp_m=temp_m+ML[area1];
                    if(temp_m>index_j-1){
                        temp_res=subpro(index_area,index_max_i,area1,index_j-1-temp_m+ML[area1]);
                        a_d[index_area][index_max_i]=m;
                        for(int area=0;area<Area;area++){
                            for(int i=0;i<N[area];i++){
                                if(a_d[area][i]==m&&area!=index_area&&i!=index_max_i){
                                    a_d[area][i]=-1;
                                }       
                            }   
                        }
                        res1[index_area][index_max_i]=temp_res[0];
                        
                        res[index_area][index_max_i]=index_j;
                        break;
                    }


                }
            }
        }
    }
    for(int area=0;area<Area;area++){
        for(int i=0;i<N[area];i++){
            if(res1[area][i]==0){
                //该用户被随机卸载
                srand(seed);
                int index_temp_1=rand()%(Area+1);
                seed+=1;
                if(index_temp_1==0){
                    temp_res=subpro(area,i,area,-1);
                    res1[area][i]=temp_res[0];
                    res[area][i]=0;
                }
                else{
                    temp_res=subpro(area,i,index_temp_1-1,0);
                    res1[area][i]=temp_res[0];
                    res[area][i]=accumulate(ML.begin(),ML.begin()+index_temp_1-1,1);
                
                }

            }
        }
    }
    // for(auto i:a){
    //     for(auto j:i){
    //         cout<<j<<' ';
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // Res=0;
  //  cout<<"m1"<<endl;
    for(auto i:res1){
        for(auto j:i){
          //  cout<<j<<' ';
            Res2+=v-p*log10(1+j);
            Res+=j;
        }
   // cout<<endl;
     }
    // cout<<Res<<endl;
    return res;
}

vector<vector<double>> exploitation1(int K,vector<int> temp_M,double& Res,double &Res2){
    vector<vector<double>> res(Area);
    vector<vector<double>> res1(Area);
    Res=0.0;
    Res2=0.0;
    for(int area=0;area<Area;area++){
        if(N[area]==0){
            continue;
        }
        res[area].resize(N[area]);
        res1[area].resize(N[area]);
        for(int i=0;i<N[area];i++){
            auto temp_res=subpro(area,i,a1_d[area][i][0],a1_d[area][i][1]);
            res1[area][i]=temp_res[0];
            int temp=a_d[area][i];
            int temp_s=0;
            int index_j=0;
            for(int k=0;k<K;k++){
                temp_s+=temp_M[k];
                if(temp_s>temp){
                    index_j=k;
                    res[area][i]=index_j;
                    break;
                }
            }
            
        }
    } 
    //cout<<"e21"<<endl;
    for(auto i:res1){
        for(auto j:i){
            //cout<<j<<' ';
            Res2+=v-p*log10(1+j);
            Res+=j;
        }
    //cout<<endl;
     }
    // cout<<Res<<endl;
    return res;
}


void experiment1(){
    vector<pair<double,double>> f_local;
    vector<pair<double,double>> p_all;
    M.clear();
    for(int i=1;i<=9;i++){
        f_local.push_back({i*0.1*pow(10,9),(1+i*0.1)*pow(10,9)});
        p_all.push_back({i*0.001,i*0.01});
    }
    int Size=7;
    int NT=10;
    vector<double> temp{5,4,2,3,2,3,4};
    int r=temp.size();
    Area=r;
    initialization_mec();
    vector<vector<double>> frequency(Size,vector<double>(r,0));

    for(int i=1;i<Size;i++){
        for(int j=0;j<r;j++){
            frequency[i][j]=temp[j]+temp[j]*(i-1);
        }
    }
    
    aa=((k1+k2)/2.0)*pow(((0.1+6)*pow(10,9)/2.0),3.0);
    cost_u=0.0;
    for(int i=0;i<Area;i++){
       cost_u=cost_u+accumulate(cost[i].begin(), cost[i].end(), 0.0)/(ML[i]);
    }
    cost_u=cost_u/Area;
    //cout<<cost_u<<endl;
    //统计期望和上下界
    int M_s=accumulate(M.begin(), M.end(), 0.0);
    //vector<vector<double>> mean_my(Area,vector<double>(M_s+Area+1));
    //double mean_max=0;
    //double mean_min=D_MAX;
    //double mean_=0;
    vector<vector<int>> T(Size-1,vector<int>(NT));
    vector<vector<double>> res_offline(Size-1);
    vector<vector<double>> res_online(Size-1);
    vector<vector<double>> time_offline(Size-1);
    vector<vector<double>> time_online(Size-1);
    vector<vector<double>> res_ucb(Size-1);
    vector<vector<double>> res_exp3(Size-1);
    vector<vector<double>> time_ucb(Size-1);
    vector<vector<double>> time_exp3(Size-1);
    vector<vector<double>> res_exp3_ix(Size-1);
    vector<vector<double>> time_exp3_ix(Size-1);
    vector<vector<double>> res_j(Size-1);
    vector<vector<double>> time_j(Size-1);
    vector<vector<double>> res_pso(Size-1);
    vector<vector<double>> time_pso(Size-1);
    vector<vector<double>> res_D(Size-1);
    vector<vector<double>> time_D(Size-1);
    // 
    double O=D_MAX;
    struct timeval sTime, eTime;
    long exeTime;
    for(int i=0;i<Size-1;i++){//每个规模
        vector<double> Nmax(Area);
        int SNmax=0;
        int Mmax=0;
        int Mmin=100;
        int SN=0;
        int SM1=0;
        int SM2=0;
        int SM1_=0;
        int SM2_=0;
        S2.clear();
        S1.clear();
        C_V.clear();
        S.clear();
        V.clear();
        V_d.clear();
        S_d.clear();
        S_d.resize(Area);
        S.resize(Area);
        C_V.resize(Area);
        V.resize(Area);
        V_d.resize(Area);
        N.clear();
        F_MAX.clear();
        S1.resize(Area);
        S2.resize(Area);
        N.resize(Area);
        F_MAX.resize(Area);
        //随机生成初始化服务器个数和用户数量
        auto f_max_temp=uniform1_int(0,10,M_s);
       // Print(f_max_temp);
        for(int j=0,jj=0;j<Area;j++){
            for(int jjj=0;jjj<M[j];jjj++){
                //cout<<f_max_temp[jj]<<' ';
                F_MAX[j].push_back(f_max_temp[jj++]);
                
            }
            //cout<<endl;
        }
        for(int j=0;j<Area;j++){
            auto temp=uniform1_int(frequency[0][j],frequency[i+1][j]-1,1);
            N[j]=temp[0];
        }
        Print(N);
        int Nm=0;
        double t2=1;
        for(int t=0;t<NT;t++){//每个阶段
            SN=0;
            SNmax=0;
            for(int j=0;j<Area;j++){
                Nmax[j]=max(Nmax[j],N[j]);
                if(Nm<Nmax[j]){
                    Nm=Nmax[j];
                }
                SNmax+=Nmax[j];
                SN+=N[j];
            }

            vector<int> temp_M;
            vector<int> temp_M_;
            vector<int> temp_M1;
            vector<int> temp_M1_;
            int min_temp=100;
            auto temp_m=uniform1(0,1,Area+1);
            auto temp_m_=uniform1(0,1,Area+1);
            double temp_m_s=accumulate(temp_m.begin(),temp_m.end(),0.0);
            double temp_m_s_=accumulate(temp_m_.begin(),temp_m_.end(),0.0);
            for(int j=0;j<Area+1;j++){
                temp_m[j]=ceil(SN*temp_m[j]/temp_m_s);
                temp_m_[j]=ceil(SNmax*temp_m_[j]/temp_m_s_);
            }
            temp_M1.push_back(SNmax);
            temp_M.push_back(SN);
            for(int j=0;j<Area;j++){
                temp_M1.push_back(SNmax);
                temp_M.push_back(SN);
                if(min_temp>SNmax){
                    min_temp=SNmax;
                }
                for(int jj=0;jj<M[j];jj++){
                    temp_M.push_back(F_MAX[j][jj]);
                    temp_M1.push_back(F_MAX[j][jj]);
                    if(min_temp>F_MAX[j][jj]){
                        min_temp=F_MAX[j][jj];
                    }
                }
                
            }
            temp_M1_.push_back(temp_m_[0]);
            temp_M_.push_back(temp_m[0]);
            for(int j=0;j<Area;j++){
                temp_M1_.push_back(temp_m_[j+1]);
                temp_M_.push_back(temp_m[j+1]);
                for(int jj=0;jj<M[j];jj++){
                    temp_M_.push_back(F_MAX[j][jj]);
                    temp_M1_.push_back(F_MAX[j][jj]);
                }
                
            }
            SM1_=accumulate(temp_M1_.begin(),temp_M1_.end(),0.0);
            SM2_=accumulate(temp_M_.begin(),temp_M_.end(),0.0);
            SM1=accumulate(temp_M1.begin(),temp_M1.end(),0.0);
            SM2=accumulate(temp_M.begin(),temp_M.end(),0.0);
            Mmax=max(Mmax,SM1);
            Mmin=min(Mmin,min_temp);
            int K=temp_M.size();

            //gap
            // vector<vector<double>> G_gap;
            // ifstream infile("mean.csv");
            // string line;
            // for (int j=0;getline(infile, line);j++)
            // {
            //     istringstream sin(line);
            //     string field;
            //     vector<double> temp;
            //     for(int jj=0;getline(sin, field, ',');jj++){
            //         temp.push_back(stod(field.c_str()));
            //     }
            //     for(int jj=0;jj<N[j];jj++){
            //         G_gap.push_back(temp);
            //     }
                
            // }
            // vector<vector<double>> G_gap_1(G_gap.size());
            // for(int j=0;j<G_gap.size();j++){
            //     for(int jj=0;jj<G_gap[0].size();jj++){
            //         for(int k=0;k<temp_M[jj];k++){
            //             G_gap_1[j].push_back(G_gap[j][jj]);
            //         }
            //     }

            // }
            // //最大值
            // lenm=G_gap_1[0].size();
            // lenn=G_gap_1.size();
            // for (int j = 0; j < lenn; ++j){
            //     for (int jj = 0; jj < lenm; ++jj){
            //         E_val[j][jj]=G_gap_1[j][jj];     
            //     }
                
            // }
            // for (int j = lenn; j < lenm; ++j){
            //     for (int jj = 0; jj < lenm; ++jj){
            //         E_val[j][jj]=0;
            //     }
                
            // }
            // double Opt;
            // auto opt_s=KM(Opt);
            // for (int j = 0; j < lenn; ++j){
            //     for (int jj = 0; jj < lenm; ++jj){
            //         if(opt_s[j][jj]==1){
            //             E_val[j][jj]=E_val[j][jj]-pow(10,-10);    
            //         } 
            //     }
                
            // }
            //double subOpt;
            //auto opt_sub=KM(subOpt);
            O=2.12453e-1;
            //cout<<O<<endl;
            double upr=3.08266e-1;
            double lowr=3.04524e-1;
            double o=2.62709e-1;
            double e1=pow(10,-20);
            double e=min(O/(4.0*SN),o/K-0.8*O/(SN*K)-e1);
            //cout<<75.0*SNmax*SNmax*Mmax*(upr-lowr)*(upr-lowr)/(8.0*O*O*Mmin)<<' '<<(exp(1)*Mmax*exp(1)*Mmax)/((3.0-exp(1))*Mmin*(3.0-exp(1))*Mmin)<<endl;
            //long T1=max(ceil(75.0*SNmax*SNmax*Mmax*(upr-lowr)*(upr-lowr)/(8.0*O*O*Mmin)),ceil((exp(1)*Mmax*exp(1)*Mmax)/((3.0-exp(1))*Mmin*(3.0-exp(1))*Mmin)));
            //long T2=ceil(SN* SM+SN* SM*upr/e);
            //cout<<e<<' '<<T1<<' '<<T2<<endl;
            //e=0.001;
            long T1=max(100.0,pow(2,t));
            long T2=10;
            T[i][t]=T1+T2+pow(2,t+1);
            for(int j=0;j<Area;j++){
                int n=S[j].size();
                int n1=S_d[j].size();
                if(n1<N[j]){
                    for(int jj=0;jj<N[j]-n;jj++){
                        S_d[j].push_back(vector<double>(K,0));
                        V_d[j].push_back(vector<double>(K,0));
                    }
                }
                if(n<Nmax[j]){
                    for(int jj=0;jj<Nmax[j]-n;jj++){
                        S[j].push_back(vector<double>(K,0));
                        V[j].push_back(vector<double>(K,0));
                        C_V[j].push_back(vector<vector<double>>(K,vector<double>(2,0)));
                        S1[j].push_back(vector<double>(K,0));
                        S2[j].push_back(vector<double>(K,0));
                    }
                }

            }
            bool flagm=true;
            bool flage=true;
            bool flagm_d=true;
            bool flage_d=true;
            a1.clear();
            R1.clear();
            B.clear();
            a1.resize(Area);
            a.clear();
            a.resize(Area);
            R1.resize(Area);
            B.resize(Area);
            a1_d.clear();
            R1_d.clear();
            B_d.clear();
            a1_d.resize(Area);
            a_d.clear();
            a_d.resize(Area);
            R1_d.resize(Area);
            B_d.resize(Area);
            int T_l=T[i][t];
            double Res=0.0;
            for(int t1=0;t1<T_l;t1++,t2++){
                //初始化数据
                initialization_user(f_local,p_all,9,9);
                //精确
                //PSO
                gettimeofday(&sTime, NULL);
                //auto res_pso1=PSO();
                //auto myG1=res_pso1[0];
                PSO_thread();
                auto myG1=res111[0];
                // for(auto j:myG1){
                //     for(auto jj:j){
                //         for(auto jjj:jj){
                //             if(jjj==D_MAX){
                //                 cout<<"error"<<endl;
                //             }
                //         }  
                //     }
                // }
                vector<vector<double>> G1_1;
                for(int j=0;j<Area;j++){
                    for(int jj=0;jj<N[j];jj++){
                        vector<double> temp_11;
                        for(int jjj=0;jjj<M_s+Area+1;jjj++){
                            for(int k=0;k<temp_M[jjj];k++){
                                temp_11.push_back(myG1[j][jj][jjj]);
                            }
                        }
                        G1_1.push_back(temp_11);
                    }
                    
                }
                

                //最大值
                lenm=G1_1[0].size();
                lenn=G1_1.size();
                for (int j = 0; j < lenn; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=10000-G1_1[j][jj]; 
                        if(E_val[j][jj]<0){
                            cout<<"error"<<endl;
                        }    
                    }
                    
                }
                for (int j = lenn; j < lenm; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=0;
                    }
                    
                }
                double Opt1=0;
                auto opt1_s=KM(Opt1);
                res_pso[i].push_back(10000*lenn-Opt1);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_pso[i].push_back(exeTime*pow(10,-3));
                
                //J-UACO
                gettimeofday(&sTime, NULL);
                double res_jj=JUACO(temp_M);
               // cout<<Res<<endl;
                res_j[i].push_back(res_jj);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_j[i].push_back(exeTime*pow(10,-3));
                //offline
                /* 测量一个事件持续的时间*/
                gettimeofday(&sTime, NULL);
                //auto res=solution_1();
                solution_11();
                auto myG=res111[0];
               // cout<<"offline"<<endl;
            //     for(auto j:myG){
            //         for(auto jj:j){
            //             for(auto jjj:jj){
            //                 // cout<<jjj<<' ';
            //                 if(jjj==D_MAX){
            //                     cout<<"error"<<endl;
            //                 }
            //             }
            //            // cout<<endl;
            //         }     
            //    }
                vector<vector<double>> G_1;
                for(int j=0;j<Area;j++){
                    for(int jj=0;jj<N[j];jj++){
                        vector<double> temp_11;
                        for(int jjj=0;jjj<M_s+Area+1;jjj++){
                            for(int k=0;k<temp_M[jjj];k++){
                                temp_11.push_back(myG[j][jj][jjj]);
                            }
                        }
                        G_1.push_back(temp_11);
                    }
                    
                }
                //最大值
                lenm=G_1[0].size();
                lenn=G_1.size();
                for (int j = 0; j < lenn; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=10000-G_1[j][jj]; 
                        if(E_val[j][jj]<0){
                            cout<<"error"<<endl;
                        }   
                        //cout<<E_val[j][jj]<<' '; 
                    }
                    //cout<<endl;
                    
                }
                for (int j = lenn; j < lenm; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=0;
                        //cout<<E_val[j][jj]<<' '; 
                    }
                    //cout<<endl; 
                    
                }
                double Opt=0;
                auto opt_s=KM(Opt);
                // for(int r=0;r<lenn;r++){
                //     for(int rr=0;rr<lenm;rr++){
                //         if(opt_s[r][rr]==1){
                //             cout<<G_1[r][rr]<<endl;
                //         }
                //     }
                // }
                //cout<<Opt<<endl;
                res_offline[i].push_back(10000*lenn-Opt);
                //cout<<1000*lenn-Opt<<endl;
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_offline[i].push_back(exeTime*pow(10,-3));
                // for(int j=0;j<Area;j++){
                //     for(int jj=0;jj<N[j];jj++){
                //         for(int jjj=0;jjj<M_s+Area+1;jjj++){
                //             double temp=v-p*log10(myG[j][jj][jjj]+1);
                //             mean_my[j][jjj]+=temp/N[j];
                //             if(mean_max<temp){
                //                 mean_max=temp;
                //             }
                //             if(mean_min>temp){
                //                 mean_min=temp;
                //             }
                //         }
                        
                //     }
                // }
                 //MUCB
                gettimeofday(&sTime, NULL);
                double Res_ucb_online=0.0;   
                double Res2_ucb_reward=0.0;  
                auto res_ucb_1=MUCB(t2,Nmax,K,temp_M1,Res_ucb_online,Res2_ucb_reward);
               // cout<<Res_ucb_online<<endl;
                res_ucb[i].push_back(Res_ucb_online);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_ucb[i].push_back(exeTime*pow(10,-3));
                //MEXP3
                gettimeofday(&sTime, NULL);
                double Res_exp_online=0.0;   
                double Res2_exp_reward=0.0;  
                auto res_exp_1=MEXP3(t2,Nmax,Nm,K,temp_M1,Res_exp_online,Res2_exp_reward);
                res_exp3[i].push_back(Res_exp_online);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_exp3[i].push_back(exeTime*pow(10,-3));
                //online  
                vector<vector<double>> res_online_temp;
                gettimeofday(&sTime, NULL);
                double Res_online=0;   
                double Res_reward=0;     
                if(t1<T1){
                   res_online_temp=exploration(Nmax,K,temp_M1_,SM1_,Nm,Res_online,Res_reward);
                }
                else if(t1>=T1&&t1<T1+T2){
                    if(flagm){
                        for(int j=0;j<Area;j++){
                            for(int jj=0;jj<N[j];jj++){
                                B[j].push_back(vector<double>(SM2,0.0));
                                a[j].push_back(-1);
                                R1[j].push_back(vector<double>(SM2,0.0));
                                for(int jjj=0,k=0;jjj<K;jjj++){
                                        for(int jjjj=0;jjjj<temp_M_[jjj];jjjj++){
                                            if(V[j][jj][jjj]!=0){
                                            R1[j][jj][k++]=S[j][jj][jjj]/V[j][jj][jjj];
                                        }
                                        else{
                                            k++;
                                        }
                                        }
                                        
                                    
                                }
                                
                            }
                        }
                            flagm=false;
                    }

                    res_online_temp=matching(e,temp_M_,SM2_,Res_online,Res_reward);
                }
               
                else if(t1>=T1+T2){                
                if(flage){
                    for(int j=0;j<Area;j++){
                        a1[j].resize(N[j],vector<double>(2,-2));
                        for(int jj=0;jj<N[j];jj++){
                            int temp=a[j][jj];
                            int temp_s=0;
                            int index_j=0;
                            for(int k=0;k<K;k++){
                                temp_s+=temp_M_[k];
                                if(temp_s>temp){
                                    index_j=k;
                                    break;
                                }
                            }
                            if(index_j==0){
                                a1[j][jj][0]=j;
                                a1[j][jj][1]=-1;
                            }   
                            else{
                                int temp_m=0;
                                for(int k=0;k<Area;k++){
                                    temp_m+=ML[k];
                                    if(temp_m>index_j-1){
                                        a1[j][jj][0]=k;
                                        a1[j][jj][1]=index_j-1-temp_m+ML[k];
                                        break;
                                    }
                                }
                            }        

                        }
                        

                    }
                flage=false;
                }
                res_online_temp=exploitation(K,temp_M_,Res_online,Res_reward);
                }
                res_online[i].push_back(Res_online);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_online[i].push_back(exeTime*pow(10,-3));

                //MEXP3_ix
                gettimeofday(&sTime, NULL);
                double Res_expix_online=0.0;   
                double Res2_expix_reward=0.0;  
                auto res_expix_1=MEXP3(t2,Nmax,Nm,K,temp_M1,Res_expix_online,Res2_expix_reward);
                res_exp3_ix[i].push_back(Res_expix_online);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_exp3_ix[i].push_back(exeTime*pow(10,-3));

                                //D-DEBO
                vector<vector<double>> res_online_d_temp;
                gettimeofday(&sTime, NULL);
                double Res_online_d=0;   
                double Res_reward_d=0;     
                if(t1<T1){
                   res_online_d_temp=exploration1(K,temp_M,SM2,Res_online_d,Res_reward_d);
                }
                else if(t1>=T1&&t1<T1+T2){
                    if(flagm_d){
                        for(int j=0;j<Area;j++){
                            for(int jj=0;jj<N[j];jj++){
                                B_d[j].push_back(vector<double>(SM2,0.0));
                                a_d[j].push_back(-1);
                                R1_d[j].push_back(vector<double>(SM2,0.0));
                                for(int jjj=0,k=0;jjj<K;jjj++){
                                        for(int jjjj=0;jjjj<temp_M[jjj];jjjj++){
                                            if(V_d[j][jj][jjj]!=0){
                                            R1_d[j][jj][k++]=S_d[j][jj][jjj]/V_d[j][jj][jjj];
                                        }
                                        else{
                                            k++;
                                        }
                                        }
                                        
                                    
                                }
                                
                            }
                        }
                            flagm_d=false;
                    }

                    res_online_d_temp=matching1(e,temp_M,SM2,Res_online_d,Res_reward_d);
                }
               
                else if(t1>=T1+T2){                
                    if(flage_d){
                    for(int j=0;j<Area;j++){
                        a1_d[j].resize(N[j],vector<double>(2,-2));
                        for(int jj=0;jj<N[j];jj++){
                            int temp=a_d[j][jj];
                            int temp_s=0;
                            int index_j=0;
                            for(int k=0;k<K;k++){
                                temp_s+=temp_M[k];
                                if(temp_s>temp){
                                    index_j=k;
                                    break;
                                }
                            }
                            if(index_j==0){
                                a1_d[j][jj][0]=j;
                                a1_d[j][jj][1]=-1;
                            }   
                            else{
                                int temp_m=0;
                                for(int k=0;k<Area;k++){
                                    temp_m+=ML[k];
                                    if(temp_m>index_j-1){
                                        a1_d[j][jj][0]=k;
                                        a1_d[j][jj][1]=index_j-1-temp_m+ML[k];
                                        break;
                                    }
                                }
                            }        

                        }
                        

                    }
                    flage_d=false;
                }
                res_online_d_temp=exploitation1(K,temp_M,Res_online_d,Res_reward_d);
                }
                res_D[i].push_back(Res_online_d);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_D[i].push_back(exeTime*pow(10,-3));
                
            }
            

            //用户数量随机生成
            for(int j=0;j<Area;j++){
                auto temp=uniform1_int(frequency[0][j],frequency[i+1][j]-1,1);
                N[j]=temp[0];
            }
            //服务器变动
           
            for(int j=0,jj=0;j<Area;j++){
                for(int jjj=0;jjj<M[j];jjj++){
                     auto temp=uniform1_int(0,10,1);
                    F_MAX[j][jjj]=temp[0];
            } 
            }
        }
    
    }
    
    // cout<<mean_min<<' '<<mean_max<<endl;
    // std::ofstream myfile1;
    // myfile1.open ("mean.csv");
    // for(int i=0;i<Area;i++){
    //     for(int j=0;j<M_s+Area+1;j++){
    //         mean_my[i][j]= mean_my[i][j]/mean_;
    //         myfile1 << mean_my[i][j];
    //         if(j<M_s+Area){
    //             myfile1<<",";
    //         }
    //     }
    //     if(i<Area-1){
    //         myfile1 << "\n";
    //     }
        
    // }

      //myfile1.close();
 
    // double mean_gap=D_MAX;
    // for(int i=0;i<Area;i++){
    //     for(int j=0;j<M_s+Area+1;j++){
    //         for(int jj=0;jj<M_s+Area+1;jj++)
    //         if(j!=jj){
    //             mean_gap=mean_gap>fabs(mean_my[i][j]-mean_my[i][jj])?fabs(mean_my[i][j]-mean_my[i][jj]):mean_gap;
    //         }
            
    //     }
    // }
     std::ofstream myfile1;
    myfile1.open ("result1_89_u.csv");
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_offline[i].begin(),res_offline[i].end(),0.0)/res_offline[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
    }
        myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_online[i].begin(),res_online[i].end(),0.0)/res_online[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
      myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_ucb[i].begin(),res_ucb[i].end(),0.0)/res_ucb[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
      myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_exp3[i].begin(),res_exp3[i].end(),0.0)/res_exp3[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
     myfile1 << "\n";
        for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_exp3_ix[i].begin(),res_exp3_ix[i].end(),0.0)/res_exp3_ix[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }

            myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_j[i].begin(),res_j[i].end(),0.0)/res_j[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
                myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_pso[i].begin(),res_pso[i].end(),0.0)/res_pso[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
                    myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(res_D[i].begin(),res_D[i].end(),0.0)/res_D[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
        myfile1 << "\n";
        for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_offline[i].begin(),time_offline[i].end(),0.0)/time_offline[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
    }
        myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_online[i].begin(),time_online[i].end(),0.0)/time_online[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
            myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_ucb[i].begin(),time_ucb[i].end(),0.0)/time_ucb[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
            myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_exp3[i].begin(),time_exp3[i].end(),0.0)/time_exp3[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
                myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_exp3_ix[i].begin(),time_exp3_ix[i].end(),0.0)/time_exp3_ix[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }

            myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_j[i].begin(),time_j[i].end(),0.0)/time_j[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
                myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_pso[i].begin(),time_pso[i].end(),0.0)/time_pso[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
                    myfile1 << "\n";
    for(int i=0;i<Size-1;i++){
            myfile1 << accumulate(time_D[i].begin(),time_D[i].end(),0.0)/time_D[i].size();
            if(i<Size-2){
                myfile1<<",";
            }
        
    }
    myfile1.close();
}

double server(int area,int i,int area1,int j,vector<vector<double>> x){
    double res=0;
    if(area==area1){
        auto p1=uniform1(p_range_all[area][i].first,p_range_all[area][i].second,1);
        if(j==0){
            res=func_MBS_SBS_111(area,i,area1,j,p1[0],f_range_mec[area][j].first,x);
        }
        else{
            auto f1=uniform1(f_range_mec[area][j].first,f_range_mec[area][j].second,1);
            res=func_MBS_SBS_111(area,i,area1,j,p1[0],f1[0],x);
        }
    }
    else{
        auto p1=uniform1(p_range_all[area][i].first,p_range_all[area][i].second,1);
        if(j==0){
            res=func_MBS_SBS_111(area,i,area1,j,p1[0],f_range_mec[area1][j].first,x);
        }
        else{
            auto f1=uniform1(f_range_mec[area1][j].first,f_range_mec[area1][j].second,1);
            res=func_MBS_SBS_111(area,i,area1,j,p1[0],f1[0],x);
        }
    }
    return res;
}

vector<vector<vector<double>>> motivtion(){
    vector<vector<vector<double>>> res(Area);
    int M_s=accumulate(M.begin(),M.end(),0);
    auto x=uniform2(pow(10,-12),5*pow(10,-9),Area,Area);
    for(int area=0;area<Area;area++){
        res[area].resize(N[area],vector<double>(M_s+Area+1));
        for(int i=0;i<N[area];i++){
            int index_g=0;  
            //
            auto f1=uniform1(f_range_local[area][i].first,f_range_local[area][i].second,1);
            res[area][i][index_g++]=func_local(area,i,f1[0]);
            for(int ii=0;ii<Area;ii++){
                //cloud
                res[area][i][index_g++]=server(area,i,ii,0,x);
                for(int k=0;k<M[ii];k++){
                    //edge
                    res[area][i][index_g++]=server(area,i,ii,k+1,x);

                }
            }
        }

    }
    return res;
}

void experiment0(){
    //用户的最优解分布
    vector<pair<double,double>> f_local;
    vector<pair<double,double>> p_all;
    int num=100;
    for(int i=1;i<=9;i++){
        f_local.push_back({i*0.1*pow(10,9),(1+i*0.1)*pow(10,9)});
        p_all.push_back({i*0.001,i*0.01});
    }
    vector<double> temp(7,num);
    N=temp;
    int r=temp.size();
    Area=r;
    initialization_mec();
    aa=((k1+k2)/2.0)*pow(((0.1+6)*pow(10,9)/2.0),3.0);
    cost_u=0.0;
    for(int i=0;i<Area;i++){
       cost_u=cost_u+accumulate(cost[i].begin(), cost[i].end(), 0.0)/(ML[i]);
    }
    cost_u=cost_u/Area;
    F_MAX.clear();
    F_MAX.resize(Area);
    int M_s=accumulate(M.begin(), M.end(), 0.0);
    auto f_max_temp=uniform1_int(0,10,M_s);
    for(int j=0,jj=0;j<Area;j++){
        for(int jjj=0;jjj<M[j];jjj++){
            F_MAX[j].push_back(f_max_temp[jj++]);     
        }
    }
    vector<vector<double>> res12(Area,vector<double>(Area+1));
    vector<vector<double>> res1(Area,vector<double>(Area+1));
    initialization_user(f_local,p_all,9,9);     
    int SN=accumulate(N.begin(),N.end(),0);
    vector<int> temp_M;  
    temp_M.push_back(SN);
    for(int j=0;j<Area;j++){
        temp_M.push_back(SN);
        for(int jj=0;jj<M[j];jj++){
            temp_M.push_back(F_MAX[j][jj]);
        }        
    }
    //
    solution_11();
    auto myG=res111[0];
                for(auto j:myG){
                    for(auto jj:j){
                        for(auto jjj:jj){
                            if(jjj==D_MAX){
                                cout<<"error"<<endl;
                            }
                        }
                    }     
               }
    auto myG1=motivtion();
    vector<vector<double>> G1_1;
    vector<vector<double>> G_1;
    for(int j=0;j<Area;j++){
        for(int jj=0;jj<N[j];jj++){
            vector<double> temp_11;
            vector<double> temp1_11;
            for(int jjj=0;jjj<M_s+Area+1;jjj++){
               // cout<<myG1[j][jj][jjj]<<' ';
                for(int k=0;k<temp_M[jjj];k++){
                    temp_11.push_back(myG[j][jj][jjj]);
                    temp1_11.push_back(myG1[j][jj][jjj]);
                }
            }
          //  cout<<endl;
            G_1.push_back(temp_11);
            G1_1.push_back(temp1_11);
        }
                    
    }    
    lenm=G_1[0].size();
    lenn=G_1.size();
    for (int j = 0; j < lenn; ++j){
        for (int jj = 0; jj < lenm; ++jj){
            E_val[j][jj]=5000-G_1[j][jj]; 
            if(E_val[j][jj]<0){
               cout<<"error1"<<endl;
            }   
                       
        }
                    
    }
    for (int j = lenn; j < lenm; ++j){
        for (int jj = 0; jj < lenm; ++jj){
            E_val[j][jj]=0;
        }         
    }
    double Opt=0;
    auto opt_s=KM(Opt);
    lenm=G1_1[0].size();
    lenn=G1_1.size();
    for (int j = 0; j < lenn; ++j){
        for (int jj = 0; jj < lenm; ++jj){
            E_val[j][jj]=5000-G1_1[j][jj]; 
            if(E_val[j][jj]<0){
                cout<<"error2"<<endl;
            }   
                       
        }
                    
    }
    for (int j = lenn; j < lenm; ++j){
        for (int jj = 0; jj < lenm; ++jj){
            E_val[j][jj]=0;
        }         
    }
    double Opt1=0;
    auto opt1_s=KM(Opt1);

    for(int j=0,area=0;j<opt_s.size();j++){
        if(j==num*(area+1)){
            area++;
        }
        int index_jj=0;
        for(int jj=0;jj<opt_s[0].size();jj++){
            if(opt_s[j][jj]==1){
                index_jj=jj;  
                //cout<<index_jj<<' ';
                break;
            }
        }
        double temp=0;
        for(int jjj=0;jjj<M_s+Area+1;jjj++){ 
            temp+=temp_M[jjj];
            if(index_jj<temp){
                index_jj=jjj;
                //cout<<index_jj<<' ';
                break;
            }                  
        }
        if(index_jj==0){
            res12[area][0]+=1;
        }
        else{
             int temp_m=0;
            for(int area1=0;area1<Area;area1++){
                temp_m+=ML[area1];
                if(temp_m>index_jj-1){
                    //cout<<area1+1;
                    res12[area][area1+1]+=1;
                    break;
                }
            }
        }
       // cout<<endl;

    }

    for(int j=0,area=0;j<opt1_s.size();j++){
        if(j==num*(area+1)){
            area++;
        }
        int index_jj=0;
        for(int jj=0;jj<opt1_s[0].size();jj++){
            if(opt1_s[j][jj]==1){
                index_jj=jj;  
                //cout<<index_jj<<' ';
                break;
            }
        }
        double temp=0;
        for(int jjj=0;jjj<M_s+Area+1;jjj++){ 
            temp+=temp_M[jjj];
            if(index_jj<temp){
                index_jj=jjj;
                //cout<<index_jj<<' ';
                break;
            }                  
        }
        if(index_jj==0){
            res1[area][0]+=1;
        }
        else{
             int temp_m=0;
            for(int area1=0;area1<Area;area1++){
                temp_m+=ML[area1];
                if(temp_m>index_jj-1){
                    //cout<<area1+1;
                    res1[area][area1+1]+=1;
                    break;
                }
            }
        }
       // cout<<endl;

    }
    std::ofstream myfile12;
    myfile12.open ("result0_1.csv");
    for(int i=0;i<res12.size();i++){
        for(int j=0;j<res12[i].size();j++){
            myfile12 << res12[i][j]/num;
            if(j<res12[i].size()-1){
                myfile12<<",";
            }
        }
        if(i<res12.size()-1) 
        myfile12 << "\n";
    }
    myfile12.close();
    std::ofstream myfile1;
    myfile1.open ("result1_1.csv");
    for(int i=0;i<res1.size();i++){
        for(int j=0;j<res1[i].size();j++){
            myfile1 << res1[i][j]/num;
            if(j<res1[i].size()-1){
                myfile1<<",";
            }
        }
        if(i<res1.size()-1) 
        myfile1 << "\n";
    }
    myfile1.close();

    





}

void experiment2(){
    vector<pair<double,double>> f_local;
    vector<pair<double,double>> p_all;
    struct timeval sTime, eTime;
    long exeTime;
    for(int i=1;i<=9;i++){
        f_local.push_back({i*0.1*pow(10,9),(1+i*0.1)*pow(10,9)});
        p_all.push_back({i*0.001,i*0.01});
    }
    int NT=10;
    vector<double> temp{5,4,2,3,2,3,4};//泊松分布
    int r=temp.size();
    Area=r;
    initialization_mec();
    aa=((k1+k2)/2.0)*pow(((0.1+6)*pow(10,9)/2.0),3.0);
    cost_u=0.0;
    for(int i=0;i<Area;i++){
       cost_u=cost_u+accumulate(cost[i].begin(), cost[i].end(), 0.0)/(ML[i]);
    }
    cost_u=cost_u/Area;
    int M_s=accumulate(M.begin(), M.end(), 0.0);
    vector<vector<double>> mean_my(Area,vector<double>(M_s+Area+1));
    double mean_=0;
    vector<int> T(NT);
    //累计与平均回报
    vector<double> Reward;
    //累计和平均遗憾
    vector<double> Regret;

    //比例
    vector<double> reward_ratio;

    vector<double> Reward_ucb;

    //累计和平均遗憾
    vector<double> Regret_ucb;

    //比例
    vector<double> reward_ratio_ucb;

    vector<double> Reward_exp;

    //累计和平均遗憾
    vector<double> Regret_exp;

    //比例
    vector<double> reward_ratio_exp;

    vector<double> Reward_exp_ix;
  
    //累计和平均遗憾
    vector<double> Regret_exp_ix;
    //比例
    vector<double> reward_ratio_exp_ix;

    vector<double> Reward_d;
    //累计和平均遗憾
    vector<double> Regret_d;

    //比例
    vector<double> reward_ratio_d;

    vector<double> Nmax(Area);
    int SNmax=0;
    int Mmax=0;
    int Mmin=100;
    int SN=0;
    int SM1=0;
    int SM2=0;
    int SM1_=0;
    int SM2_=0;
    S_d.clear();
    S_d.resize(Area);
    V_d.clear();
    V_d.resize(Area);
    S.clear();
    C_V.clear();
    S1.clear();
    S2.clear();
    V.clear();
    S.resize(Area);
    S2.resize(Area);
    S1.resize(Area);
    V.resize(Area);
    C_V.resize(Area);
    N.clear();
    F_MAX.clear();
    N.resize(Area);
    F_MAX.resize(Area);
    vector<vector<double>> N_T(NT,vector<double>(Area));
    //随机生成初始化服务器个数和用户数量
    for(int j=0;j<Area;j++){
        auto temp_u=poisson(temp[j],NT);
        for(int jj=0;jj<NT;jj++){
            N_T[jj][j]=temp_u[jj];
        }
        
    }
    auto f_max_temp=uniform1_int(0,10,M_s);
    for(int j=0,jj=0;j<Area;j++){
        for(int jjj=0;jjj<M[j];jjj++)
        F_MAX[j].push_back(f_max_temp[jj++]);
    }
    int Nm=0;
    double t2=1;
    vector<vector<double>> res_offline(NT);
    vector<vector<double>> res_online(NT);
    vector<vector<double>> time_offline(NT);
    vector<vector<double>> time_online(NT);
    vector<vector<double>> res_ucb(NT);
    vector<vector<double>> res_exp3(NT);
    vector<vector<double>> time_ucb(NT);
    vector<vector<double>> time_exp3(NT);
    vector<vector<double>> res_exp3_ix(NT);
    vector<vector<double>> time_exp3_ix(NT);
    vector<vector<double>> res_j(NT);
    vector<vector<double>> time_j(NT);
    vector<vector<double>> res_pso(NT);
    vector<vector<double>> time_pso(NT);
    vector<vector<double>> res_d(NT);
    vector<vector<double>> time_d(NT);
    for(int t=0;t<NT;t++){//每个阶段
        N=N_T[t];
        SN=0;
        SNmax=0;
        for(int j=0;j<Area;j++){
            Nmax[j]=max(Nmax[j],N[j]);
            if(Nm<Nmax[j]){
                Nm=Nmax[j];
            }
            SNmax+=Nmax[j];
            SN+=N[j];
        }

        vector<int> temp_M;//探索的范围，不是实际上限
        vector<int> temp_M1; 
        vector<int> temp_M_;//探索的范围，不是实际上限
        vector<int> temp_M1_; 
        int min_temp=100;    
        temp_M1.push_back(SNmax);
        temp_M.push_back(SN);
        for(int j=0;j<Area;j++){
            temp_M.push_back(SN);
            temp_M1.push_back(SNmax);
            if(min_temp>SNmax){
                min_temp=SNmax;
            }
            for(int jj=0;jj<M[j];jj++){
                temp_M.push_back(F_MAX[j][jj]);
                temp_M1.push_back(F_MAX[j][jj]);
                if(min_temp>F_MAX[j][jj]){
                    min_temp=F_MAX[j][jj];
                }
            }
                
        }
        auto temp_m=uniform1(0,1,Area+1);
        auto temp_m_=uniform1(0,1,Area+1);
        double temp_m_s=accumulate(temp_m.begin(),temp_m.end(),0.0);
        double temp_m_s_=accumulate(temp_m_.begin(),temp_m_.end(),0.0);
        for(int j=0;j<Area+1;j++){
            temp_m[j]=ceil(SN*temp_m[j]/temp_m_s);
            temp_m_[j]=ceil(SNmax*temp_m_[j]/temp_m_s_);
        }
        temp_M1_.push_back(temp_m_[0]);
        temp_M_.push_back(temp_m[0]);
        for(int j=0;j<Area;j++){
            temp_M1_.push_back(temp_m_[j+1]);
            temp_M_.push_back(temp_m[j+1]);
            for(int jj=0;jj<M[j];jj++){
                temp_M_.push_back(F_MAX[j][jj]);
                temp_M1_.push_back(F_MAX[j][jj]);
            }
                
        }
        SM1=accumulate(temp_M1.begin(),temp_M1.end(),0.0);
        SM2=accumulate(temp_M.begin(),temp_M.end(),0.0);
        SM1_=accumulate(temp_M1_.begin(),temp_M1_.end(),0.0);
        SM2_=accumulate(temp_M_.begin(),temp_M_.end(),0.0);
        Mmax=max(Mmax,SM1);
        Mmin=min(Mmin,min_temp);
        int K=temp_M.size();
        double O=2.12453e-1;
        double upr=3.08266e-1;
        double lowr=3.04524e-1;
        double o=2.62709e-1;
        double e1=pow(10,-20);
        double e=min(O/(4.0*SN),o/K-0.8*O/(SN*K)-e1);
        long T1=max(50.0,pow(2,t));
        long T2=10;
        T[t]=T1+T2+pow(2,t+1);
        mean_+=T[t];
        for(int j=0;j<Area;j++){
            int n=S[j].size();
            int n1=S_d[j].size();
            if(n1<N[j]){
                for(int jj=0;jj<N[j]-n;jj++){
                    S_d[j].push_back(vector<double>(K,0));
                    V_d[j].push_back(vector<double>(K,0));
                }
            }
            if(n<Nmax[j]){
                for(int jj=0;jj<Nmax[j]-n;jj++){
                    S[j].push_back(vector<double>(K,0));
                    V[j].push_back(vector<double>(K,0));
                    C_V[j].push_back(vector<vector<double>>(K,vector<double>(2,0)));
                    S1[j].push_back(vector<double>(K,0));
                    S2[j].push_back(vector<double>(K,0));
                }
            }

        }
        //此阶段的最优回报
        vector<vector<double>> G_gap;//每种决策的回报
        ifstream infile("mean.csv");
        string line;
        for (int j=0;getline(infile, line);j++)
        {
            istringstream sin(line);
            string field;
            vector<double> temp;
            for(int jj=0;getline(sin, field, ',');jj++){
                temp.push_back(stod(field.c_str()));
                if(stod(field.c_str())>1){
                    cout<<"error1"<<endl;
                }
            } 
                G_gap.push_back(temp);
            
        }
        vector<vector<double>> G_gap_0;
        for(int j=0;j<Area;j++){
            for(int jj=0;jj<N[j];jj++){
                G_gap_0.push_back(G_gap[j]);
            }
        }
        vector<vector<double>> G_gap_1(G_gap_0.size());
        for(int j=0;j<G_gap_0.size();j++){
            for(int jj=0;jj<G_gap_0[0].size();jj++){
                for(int k=0;k<temp_M[jj];k++){
                    G_gap_1[j].push_back(G_gap_0[j][jj]);
                }
            }

        }
            //最大值
        lenm=G_gap_1[0].size();
        lenn=G_gap_1.size();
        for (int j = 0; j < lenn; ++j){
            for (int jj = 0; jj < lenm; ++jj){
                E_val[j][jj]=G_gap_1[j][jj];     
            }
                
        }
        for (int j = lenn; j < lenm; ++j){
            for (int jj = 0; jj < lenm; ++jj){
                E_val[j][jj]=0;
            }
                
        }
        double Opt;
        auto opt_s=KM(Opt);
        bool flagm=true;
        bool flage=true;
        bool flagm_d=true;
        bool flage_d=true;
        a1.clear();
        R1.clear();
        B.clear();
        a1.resize(Area);
        a.clear();
        a.resize(Area);
        R1.resize(Area);
        B.resize(Area);
        a1_d.clear();
        R1_d.clear();
        B_d.clear();
        a1_d.resize(Area);
        a_d.clear();
        a_d.resize(Area);
        R1_d.resize(Area);
        B_d.resize(Area);
        int T_l=T[t];
        double Res=0.0;
        for(int t1=0;t1<T_l;t1++,t2++){
            //初始化数据
            initialization_user(f_local,p_all,9,9);   
            //pso
            gettimeofday(&sTime, NULL);
            PSO_thread();
            auto myG1=res111[0];
                // for(auto j:myG1){
                //     for(auto jj:j){
                //         for(auto jjj:jj){
                //             if(jjj==D_MAX){
                //                 cout<<"error"<<endl;
                //             }

                //         }
                        
                //     }

                    
                // }
                vector<vector<double>> G1_1;
                for(int j=0;j<Area;j++){
                    for(int jj=0;jj<N[j];jj++){
                        vector<double> temp_11;
                        for(int jjj=0;jjj<M_s+Area+1;jjj++){
                            for(int k=0;k<temp_M[jjj];k++){
                                temp_11.push_back(myG1[j][jj][jjj]);
                            }
                        }
                        G1_1.push_back(temp_11);
                    }
                    
                }
                //最大值
                lenm=G1_1[0].size();
                lenn=G1_1.size();
                for (int j = 0; j < lenn; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=10000-G1_1[j][jj]; 
                        if(E_val[j][jj]<0){
                            cout<<"error1"<<endl;
                        }    
                    }
                    
                }
                for (int j = lenn; j < lenm; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=0;
                    }
                    
                }
                double Opt11=0;
                auto opt1_s=KM(Opt11);
                res_pso[t].push_back(10000*lenn-Opt11);
            gettimeofday(&eTime, NULL);
            exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
            time_pso[t].push_back(exeTime*pow(10,-3));
            //J
            gettimeofday(&sTime, NULL);
            double res_jj=JUACO(temp_M);
            // cout<<Res<<endl;
            res_j[t].push_back(res_jj);
            gettimeofday(&eTime, NULL);
            exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
            time_j[t].push_back(exeTime*pow(10,-3));
            gettimeofday(&sTime, NULL);
            solution_11();
            auto myG=res111[0];
                vector<vector<double>> G;
                for(auto j:myG){
                    for(auto jj:j){
                       G.push_back(jj);
                    }
                    
                }
                vector<vector<double>> G_1(G.size());
                for(int j=0;j<G.size();j++){
                    for(int jj=0;jj<G[0].size();jj++){
                        for(int k=0;k<temp_M[jj];k++){
                            G_1[j].push_back(G[j][jj]);
                        }
                    }

                }

                //最大值
                lenm=G_1[0].size();
                lenn=G_1.size();
                
                for (int j = 0; j < lenn; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=10000-G_1[j][jj]; 
                        if(E_val[j][jj]<0){
                            cout<<"error2"<<endl;
                        }    
                    }
                    
                }
                for (int j = lenn; j < lenm; ++j){
                    for (int jj = 0; jj < lenm; ++jj){
                        E_val[j][jj]=0;
                    }
                    
                }
                double Opt1=0;
                auto opt_s1=KM(Opt1);
                res_offline[t].push_back(10000*lenn-Opt1);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_offline[t].push_back(exeTime*pow(10,-3));

                
            for(int j=0;j<Area;j++){
                for(int jj=0;jj<N[j];jj++){
                    for(int jjj=0;jjj<M_s+Area+1;jjj++){
                        double temp=v-p*log10(myG[j][jj][jjj]+1);
                        mean_my[j][jjj]+=temp/N[j];
                    }
                    
                }
            }

            // vector<vector<double>> G;
            //     for(auto j:myG){
            //         for(auto jj:j){
            //            G.push_back(jj);
            //         }
                    
            //     }
            //     vector<vector<double>> G_1(G.size());
            //     //Print(N);
            //     for(int j=0;j<G.size();j++){
            //         for(int jj=0;jj<G[0].size();jj++){
            //             for(int k=0;k<temp_M[jj];k++){
            //                 G_1[j].push_back(G[j][jj]);
            //             }
            //         }
            //     }
            //     //最大值
            //     lenm=G_1[0].size();
            //     lenn=G_1.size();
            //     for (int j = 0; j < lenn; ++j){
            //         for (int jj = 0; jj < lenm; ++jj){
            //             E_val[j][jj]=v-p*log10(G_1[j][jj]+1); 
            //             if(E_val[j][jj]<0){
            //                 cout<<"error3"<<endl;
            //             }    
            //         }
                    
            //     }
            //     for (int j = lenn; j < lenm; ++j){
            //         for (int jj = 0; jj < lenm; ++jj){
            //             E_val[j][jj]=0;
            //         }
                    
            //     }
            //     double Opt_r=0;
            //     auto opt_s_r=KM(Opt_r);
            //online  
      
            gettimeofday(&sTime, NULL);
            double Res_online=0;  
            double Res_reward=0.0;
            vector<vector<double>> res_online_temp;      
            if(t1<T1){
                res_online_temp=exploration(Nmax,K,temp_M1_,SM1_,Nm,Res_online,Res_reward);
            }
            else if(t1>=T1&&t1<T1+T2){
                if(flagm){
                    for(int j=0;j<Area;j++){
                        for(int jj=0;jj<N[j];jj++){
                            B[j].push_back(vector<double>(SM2,0.0));
                            a[j].push_back(-1);
                            R1[j].push_back(vector<double>(SM2,0.0));
                            for(int jjj=0,k=0;jjj<K;jjj++){
                                    for(int jjjj=0;jjjj<temp_M_[jjj];jjjj++){
                                        if(V[j][jj][jjj]!=0){
                                        R1[j][jj][k++]=S[j][jj][jjj]/V[j][jj][jjj];
                                    }
                                    else{
                                        k++;
                                    }
                                }      
                            }
                                
                        }
                    }
                        flagm=false;
                }
                res_online_temp=matching(e,temp_M_,SM2_,Res_online,Res_reward);
            }
            else if(t1>=T1+T2){                
                if(flage){
                    for(int j=0;j<Area;j++){
                        a1[j].resize(N[j],vector<double>(2,-2));
                        for(int jj=0;jj<N[j];jj++){
                            int temp=a[j][jj];
                            int temp_s=0;
                            int index_j=0;
                            for(int k=0;k<K;k++){
                                temp_s+=temp_M_[k];
                                if(temp_s>temp){
                                    index_j=k;
                                    break;
                                }
                            }
                            if(index_j==0){
                                a1[j][jj][0]=j;
                                a1[j][jj][1]=-1;
                            }   
                            else{
                                int temp_m=0;
                                for(int k=0;k<Area;k++){
                                    temp_m+=ML[k];
                                    if(temp_m>index_j-1){
                                        a1[j][jj][0]=k;
                                        a1[j][jj][1]=index_j-1-temp_m+ML[k];
                                        break;
                                    }
                                }
                            }        

                        }
                    }
                    flage=false;
                }
                res_online_temp=exploitation(K,temp_M_,Res_online,Res_reward);
            }
            res_online[t].push_back(Res_online);
            if(Reward.empty()){
                Reward.push_back(Res_reward/SN);
            }
            else{
                Reward.push_back(Res_reward/SN+Reward.back());
            }

            reward_ratio.push_back((Reward.back()/Reward.size())*SN/Opt);
            //此时的实际期望
            double opt1=0.0;
            for(int i=0;i<Area;i++){
                for(int j=0;j<N[i];j++){
                if(G_gap[i][res_online_temp[i][j]]>1){
                    cout<<"error4"<<endl;
                }
                    opt1+=G_gap[i][res_online_temp[i][j]];
                }
            }
            //遗憾
            
            if(Regret.empty()){
                Regret.push_back(Opt/SN-opt1/SN);
            }
            else{
                Regret.push_back(Opt/SN-opt1/SN+Regret.back());
            }
            gettimeofday(&eTime, NULL);
            exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
            time_online[t].push_back(exeTime*pow(10,-3)); 
             //MUCB
            gettimeofday(&sTime, NULL);
            double Res_ucb_online=0.0;   
            double Res2_ucb_reward=0.0;  
            auto res_ucb_1=MUCB(t2,Nmax,K,temp_M1,Res_ucb_online,Res2_ucb_reward);
            res_ucb[t].push_back(Res_ucb_online);
            if(Reward_ucb.empty()){
                Reward_ucb.push_back(Res2_ucb_reward/SN);
            }
            else{
                Reward_ucb.push_back(Res2_ucb_reward/SN+Reward_ucb.back());
            }
            reward_ratio_ucb.push_back((Reward_ucb.back()/Reward_ucb.size())*SN/Opt);
            //此时的实际期望
            double opt1_ucb=0.0;
            for(int i=0;i<Area;i++){
                for(int j=0;j<N[i];j++){
                if(G_gap[i][res_ucb_1[i][j]]>1){
                    cout<<"error5"<<endl;
                }
                    opt1_ucb+=G_gap[i][res_ucb_1[i][j]];
                }
            }
            //遗憾
            
            if(Regret_ucb.empty()){
                Regret_ucb.push_back(Opt/SN-opt1_ucb/SN);
            }
            else{
                Regret_ucb.push_back(Opt/SN-opt1_ucb/SN+Regret_ucb.back());
            }
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_ucb[t].push_back(exeTime*pow(10,-3));
                //MEXP3
                gettimeofday(&sTime, NULL);
                double Res_exp_online=0.0;   
                double Res2_exp_reward=0.0;  
                auto res_exp_1=MEXP3(t2,Nmax,Nm,K,temp_M1,Res_exp_online,Res2_exp_reward);
                if(Reward_exp.empty()){
                Reward_exp.push_back(Res2_exp_reward/SN);
            }
            else{
                Reward_exp.push_back(Res2_exp_reward/SN+Reward_exp.back());
            }
            reward_ratio_exp.push_back((Reward_exp.back()/Reward_exp.size())*SN/Opt);
            //此时的实际期望
            double opt1_exp=0.0;
            for(int i=0;i<Area;i++){
                for(int j=0;j<N[i];j++){
                if(G_gap[i][res_exp_1[i][j]]>1){
                    cout<<"error3"<<endl;
                }
                    opt1_exp+=G_gap[i][res_exp_1[i][j]];
                }
            }
            //遗憾
            
            if(Regret_exp.empty()){
                Regret_exp.push_back(Opt/SN-opt1_exp/SN);
            }
            else{
                Regret_exp.push_back(Opt/SN-opt1_exp/SN+Regret_exp.back());
            }
                res_exp3[t].push_back(Res_exp_online);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_exp3[t].push_back(exeTime*pow(10,-3));


                //MEXP3_ix
                gettimeofday(&sTime, NULL);
                double Res_expix_online=0.0;   
                double Res2_expix_reward=0.0;  
                auto res_expix_1=MEXP3(t2,Nmax,Nm,K,temp_M1,Res_expix_online,Res2_expix_reward);
                if(Reward_exp_ix.empty()){
                Reward_exp_ix.push_back(Res2_expix_reward/SN);
            }
            else{
                Reward_exp_ix.push_back(Res2_expix_reward/SN+Reward_exp_ix.back());
            }
            reward_ratio_exp_ix.push_back((Reward_exp_ix.back()/Reward_exp_ix.size())*SN/Opt);
            //此时的实际期望
            double opt1_exp_ix=0.0;
            for(int i=0;i<Area;i++){
                for(int j=0;j<N[i];j++){
                if(G_gap[i][res_expix_1[i][j]]>1){
                    cout<<"error4"<<endl;
                }
                    opt1_exp_ix+=G_gap[i][res_expix_1[i][j]];
                }
            }
            //遗憾
            
            if(Regret_exp_ix.empty()){
                Regret_exp_ix.push_back(Opt/SN-opt1_exp_ix/SN);
            }
            else{
                Regret_exp_ix.push_back(Opt/SN-opt1_exp_ix/SN+Regret_exp_ix.back());
            }
                res_exp3_ix[t].push_back(Res_expix_online);
                gettimeofday(&eTime, NULL);
                exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
                time_exp3_ix[t].push_back(exeTime*pow(10,-3));

                 //DEBO
            gettimeofday(&sTime, NULL);
            double Res_online_d=0;  
            double Res_reward_d=0.0;
            vector<vector<double>> res_online_temp_d;      
            if(t1<T1){
                res_online_temp_d=exploration1(K,temp_M,SM2,Res_online_d,Res_reward_d);
            }
            else if(t1>=T1&&t1<T1+T2){
                if(flagm_d){
                    for(int j=0;j<Area;j++){
                        for(int jj=0;jj<N[j];jj++){
                            B_d[j].push_back(vector<double>(SM2,0.0));
                            a_d[j].push_back(-1);
                            R1_d[j].push_back(vector<double>(SM2,0.0));
                            for(int jjj=0,k=0;jjj<K;jjj++){
                                    for(int jjjj=0;jjjj<temp_M[jjj];jjjj++){
                                        if(V_d[j][jj][jjj]!=0){
                                        R1_d[j][jj][k++]=S_d[j][jj][jjj]/V_d[j][jj][jjj];
                                    }
                                    else{
                                        k++;
                                    }
                                }      
                            }
                                
                        }
                    }
                        flagm_d=false;
                }
                res_online_temp_d=matching(e,temp_M,SM2,Res_online_d,Res_reward_d);
            }
            else if(t1>=T1+T2){                
                if(flage_d){
                    for(int j=0;j<Area;j++){
                        a1_d[j].resize(N[j],vector<double>(2,-2));
                        for(int jj=0;jj<N[j];jj++){
                            int temp=a_d[j][jj];
                            int temp_s=0;
                            int index_j=0;
                            for(int k=0;k<K;k++){
                                temp_s+=temp_M[k];
                                if(temp_s>temp){
                                    index_j=k;
                                    break;
                                }
                            }
                            if(index_j==0){
                                a1_d[j][jj][0]=j;
                                a1_d[j][jj][1]=-1;
                            }   
                            else{
                                int temp_m=0;
                                for(int k=0;k<Area;k++){
                                    temp_m+=ML[k];
                                    if(temp_m>index_j-1){
                                        a1_d[j][jj][0]=k;
                                        a1_d[j][jj][1]=index_j-1-temp_m+ML[k];
                                        break;
                                    }
                                }
                            }        

                        }
                    }
                    flage_d=false;
                }
                res_online_temp_d=exploitation(K,temp_M,Res_online_d,Res_reward_d);
            }
            res_d[t].push_back(Res_online_d);
            if(Reward_d.empty()){
                Reward_d.push_back(Res_reward_d/SN);
            }
            else{
                Reward_d.push_back(Res_reward_d/SN+Reward_d.back());
            }

            reward_ratio_d.push_back((Reward_d.back()/Reward_d.size())*SN/Opt);
            //此时的实际期望
            double opt1_d=0.0;
            for(int i=0;i<Area;i++){
                for(int j=0;j<N[i];j++){
                if(G_gap[i][res_online_temp_d[i][j]]>1){
                    cout<<"error4"<<endl;
                }
                    opt1_d+=G_gap[i][res_online_temp_d[i][j]];
                }
            }
            //遗憾
            
            if(Regret_d.empty()){
                Regret_d.push_back(Opt/SN-opt1_d/SN);
            }
            else{
                Regret_d.push_back(Opt/SN-opt1_d/SN+Regret_d.back());
            }
            gettimeofday(&eTime, NULL);
            exeTime = (eTime.tv_sec-sTime.tv_sec)*1000000+(eTime.tv_usec-sTime.tv_usec);
            time_d[t].push_back(exeTime*pow(10,-3)); 
        }        
        //服务器变动
        
        for(int j=0,jj=0;j<Area;j++){
            for(int jjj=0;jjj<M[j];jjj++){
                auto temp=uniform1_int(0,10,1);
                F_MAX[j][jjj]=temp[0];
        } 
        }
    }

   // Print(Regret);

    std::ofstream myfile1;
    myfile1.open ("mean.csv");
    for(int i=0;i<Area;i++){
        for(int j=0;j<M_s+Area+1;j++){
            mean_my[i][j]= mean_my[i][j]/mean_;
            myfile1 << mean_my[i][j];
            if(j<M_s+Area){
                myfile1<<",";
            }
        }
        if(i<Area-1){
            myfile1 << "\n";
        }
    }
    myfile1.close();
    std::ofstream myfile;
    myfile.open ("result9_11.csv");
    for(int i=0;i<NT;i++){
            myfile << accumulate(res_offline[i].begin(),res_offline[i].end(),0.0)/res_offline[i].size();
            if(i<NT-1){
                myfile<<",";
            }
    }
        myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(res_online[i].begin(),res_online[i].end(),0.0)/res_online[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
      myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(res_ucb[i].begin(),res_ucb[i].end(),0.0)/res_ucb[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
      myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(res_exp3[i].begin(),res_exp3[i].end(),0.0)/res_exp3[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
     myfile << "\n";
        for(int i=0;i<NT;i++){
            myfile << accumulate(res_exp3_ix[i].begin(),res_exp3_ix[i].end(),0.0)/res_exp3_ix[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
            myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(res_j[i].begin(),res_j[i].end(),0.0)/res_j[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
                myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(res_pso[i].begin(),res_pso[i].end(),0.0)/res_pso[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
            myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(res_d[i].begin(),res_d[i].end(),0.0)/res_d[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
        myfile << "\n";
        for(int i=0;i<NT;i++){
            myfile << accumulate(time_offline[i].begin(),time_offline[i].end(),0.0)/time_offline[i].size();
            if(i<NT-1){
                myfile<<",";
            }
    }
        myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(time_online[i].begin(),time_online[i].end(),0.0)/time_online[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
            myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(time_ucb[i].begin(),time_ucb[i].end(),0.0)/time_ucb[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
            myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(time_exp3[i].begin(),time_exp3[i].end(),0.0)/time_exp3[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
                myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(time_exp3_ix[i].begin(),time_exp3_ix[i].end(),0.0)/time_exp3_ix[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
    
            myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(time_j[i].begin(),time_j[i].end(),0.0)/time_j[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
            myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(time_pso[i].begin(),time_pso[i].end(),0.0)/time_pso[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }
                myfile << "\n";
    for(int i=0;i<NT;i++){
            myfile << accumulate(time_d[i].begin(),time_d[i].end(),0.0)/time_d[i].size();
            if(i<NT-1){
                myfile<<",";
            }
        
    }

    myfile.close();

    myfile1.open ("result9_21.csv");
        for(int j=0;j<Regret.size();j++){
            myfile1 <<Regret[j];
            if(j<Regret.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Regret.size();j++){
            myfile1 <<Regret[j]/(j+1.0);
            if(j<Regret.size()-1){
                myfile1<<",";
            }
        }
         myfile1 << "\n";
        for(int j=0;j<Reward.size();j++){
            myfile1 <<Reward[j];
            if(j<Reward.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Reward.size();j++){
            myfile1 <<Reward[j]/(j+1.0);
            if(j<Reward.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<reward_ratio.size();j++){
            myfile1 <<reward_ratio[j];
            if(j<reward_ratio.size()-1){
                myfile1<<",";
            }
        }
        //ucb
        
         myfile1 << "\n";
        for(int j=0;j<Regret_ucb.size();j++){
            myfile1 <<Regret_ucb[j];
            if(j<Regret_ucb.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Regret_ucb.size();j++){
            myfile1 <<Regret_ucb[j]/(j+1.0);
            if(j<Regret_ucb.size()-1){
                myfile1<<",";
            }
        }
         myfile1 << "\n";
        for(int j=0;j<Reward_ucb.size();j++){
            myfile1 <<Reward_ucb[j];
            if(j<Reward_ucb.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Reward_ucb.size();j++){
            myfile1 <<Reward_ucb[j]/(j+1.0);
            if(j<Reward_ucb.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<reward_ratio_ucb.size();j++){
            myfile1 <<reward_ratio_ucb[j];
            if(j<reward_ratio_ucb.size()-1){
                myfile1<<",";
            }
        }
        //exp
        myfile1 << "\n";
        for(int j=0;j<Regret_exp.size();j++){
            myfile1 <<Regret_exp[j];
            if(j<Regret_exp.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Regret_exp.size();j++){
            myfile1 <<Regret_exp[j]/(j+1.0);
            if(j<Regret_exp.size()-1){
                myfile1<<",";
            }
        }
         myfile1 << "\n";
        for(int j=0;j<Reward_exp.size();j++){
            myfile1 <<Reward_exp[j];
            if(j<Reward_exp.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Reward_exp.size();j++){
            myfile1 <<Reward_exp[j]/(j+1.0);
            if(j<Reward_exp.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<reward_ratio_exp.size();j++){
            myfile1 <<reward_ratio_exp[j];
            if(j<reward_ratio_exp.size()-1){
                myfile1<<",";
            }
        }
        //
         myfile1 << "\n";
        for(int j=0;j<Regret_exp_ix.size();j++){
            myfile1 <<Regret_exp_ix[j];
            if(j<Regret_exp_ix.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Regret_exp_ix.size();j++){
            myfile1 <<Regret_exp_ix[j]/(j+1.0);
            if(j<Regret_exp_ix.size()-1){
                myfile1<<",";
            }
        }
         myfile1 << "\n";
        for(int j=0;j<Reward_exp_ix.size();j++){
            myfile1 <<Reward_exp_ix[j];
            if(j<Reward_exp_ix.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Reward_exp_ix.size();j++){
            myfile1 <<Reward_exp_ix[j]/(j+1.0);
            if(j<Reward_exp_ix.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<reward_ratio_exp_ix.size();j++){
            myfile1 <<reward_ratio_exp_ix[j];
            if(j<reward_ratio_exp_ix.size()-1){
                myfile1<<",";
            }
        }
        //
         myfile1 << "\n";
        for(int j=0;j<Regret_d.size();j++){
            myfile1 <<Regret_d[j];
            if(j<Regret_d.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Regret_d.size();j++){
            myfile1 <<Regret_d[j]/(j+1.0);
            if(j<Regret_d.size()-1){
                myfile1<<",";
            }
        }
         myfile1 << "\n";
        for(int j=0;j<Reward_d.size();j++){
            myfile1 <<Reward_d[j];
            if(j<Reward_d.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<Reward_d.size();j++){
            myfile1 <<Reward_d[j]/(j+1.0);
            if(j<Reward_d.size()-1){
                myfile1<<",";
            }
        }
        myfile1 << "\n";
        for(int j=0;j<reward_ratio_d.size();j++){
            myfile1 <<reward_ratio_d[j];
            if(j<reward_ratio_d.size()-1){
                myfile1<<",";
            }
        }
    myfile1.close();



}

void performance(){
    ifstream infile("result1_89_u.csv");
    vector<vector<double>> d(16);
    string line;
    for (int i=0;getline(infile, line);i++)
    {
        istringstream sin(line);
        string field;
        while (getline(sin, field, ',')) 
        {
            d[i].push_back(stod(field.c_str()));
        }
       
    }
    int n=d[0].size();
    double res=0.0,res1=0.0,res2=0.0,res3=0.0;
    for(int j=0;j<n;j++){
        res+=(d[5][j]-d[0][j])/d[5][j];
        res1+=(d[6][j]-d[0][j])/d[6][j];
        res2+=(d[13][j]-d[8][j])/d[13][j];
        res3+=(d[14][j]-d[8][j])/d[14][j];
    }
    cout<<"offline_qoe:"<<res/n<<" "<<res1/n<<endl;
    cout<<"offline_time:"<<res2/n<<" "<<res3/n<<endl;
    res=0.0;
    res1=0.0;
    res2=0.0;
    res3=0.0;
    for(int j=0;j<n;j++){
        res+=(d[2][j]-d[1][j])/d[2][j];
        res1+=(d[3][j]-d[1][j])/d[3][j];
        res2+=(d[4][j]-d[1][j])/d[4][j];
        res3+=(d[7][j]-d[1][j])/d[7][j];
    }
    cout<<"online_qoe:"<<res/n<<" "<<res1/n<<" "<<res2/n<<" "<<res3/n<<endl;

    //time
    ifstream infile1("result9_11.csv");
    vector<vector<double>> d1(16);
    string line1;
    for (int i=0;getline(infile1, line1);i++)
    {
        istringstream sin1(line1);
        string field1;
        while (getline(sin1, field1, ',')) 
        {
            d1[i].push_back(stod(field1.c_str()));
        }
       
    }
    n=d1[0].size();
    res=0.0,res1=0.0,res2=0.0,res3=0.0;
    for(int j=0;j<n;j++){
        res+=(d1[5][j]-d1[0][j])/d1[5][j];
        res1+=(d1[6][j]-d1[0][j])/d1[6][j];
        res2+=(d1[13][j]-d1[8][j])/d1[13][j];
        res3+=(d1[14][j]-d1[8][j])/d1[14][j];
    }
    cout<<"offline_qoe:"<<res/n<<" "<<res1/n<<endl;
    cout<<"offline_time:"<<res2/n<<" "<<res3/n<<endl;
    res=0.0;
    res1=0.0;
    res2=0.0;
    res3=0.0;
    for(int j=0;j<n;j++){
        res+=(d1[2][j]-d1[1][j])/d1[2][j];
        res1+=(d1[3][j]-d1[1][j])/d1[3][j];
        res2+=(d1[4][j]-d1[1][j])/d1[4][j];
        res3+=(d1[7][j]-d1[1][j])/d1[7][j];
    }
    cout<<"online_qoe:"<<res/n<<" "<<res1/n<<" "<<res2/n<<" "<<res3/n<<endl;

  
}


int main()
{

    fibonacci_creat();
    FL=fibonacci_num.size();
   // experiment0();
   // experiment1();
    experiment2();
     performance();
    return 1;
}

