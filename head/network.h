#pragma once
#include <qpOASES.hpp>
#include<string>
#include<fstream>

USING_NAMESPACE_QPOASES
using namespace std;
class master;
class worker
{
public:
	worker() {};
	worker(int m, int n, float rho_input, int *master_id_input, int worker_tag, ifstream &ReadFile);
	~worker() {};
	void updateH(float rho);
	void updateg(float rho);
	void updatex(float rho);
	void updatelambda(float rho);
	//void receivez(master *masters);
	void updateSend(float rho);
	float getR_2(master* masters);
	//float getS_2();
	float getPri();
	float getDual();
	void receiveZ(master *masters, float over_relax);

	real_t *send_val;
	real_t *z;
	real_t *x;
	real_t *x_temp;
	real_t *lambda;
	int m, n;

private:
	
	bool hot_flag;
	real_t *H_temp;
	real_t *g_temp;
	real_t *H;
	real_t *g;
	real_t *A;
	real_t *lbA;
	real_t *ubA;
	real_t *lb;
	real_t *ub;
	SQProblem model;
	//float rho;
	
	

	int worker_tag;
	int *master_id;
	
};

class master
{
public:
	master() {};
	master(int m, int worker_num_input, int *con_vec_input, int* worker_pos, ifstream &ReadFile);
	~master() {};
	void updatez(worker *workers);
	float getS_2(worker *workers);
	float getLambda(worker *workers);
	float getPri();
	void printZ();

	real_t *z;
private:
	int m;
	int worker_num;
	int *workers_id;
	int *worker_pos;
	real_t *z_temp;

};

class ADMM
{
public:
	ADMM(int M, int W, int m, int n, float rho_input, 
		int *net_vec_input, int* worker_nums, int *worker_ids, int* worker_pos, 
		float eps_rel, float eps_abs, float tau_inc, float tau_dec, float mu, int cumu, float over_relax,
		string worker_file, string master_file);
	~ADMM() {};
	void updatePri();
	void updateDual();
	bool check_conv(int tag, int max_step, int output_step);
	void updateRho();
	void solve(int max_step, int output_step, int thread_num);
	real_t * getSolution();

private:
	master *masters;
	worker *workers;
	float eps_rel, eps_abs, eps_pri, eps_dual, R_2, S_2;
	float tau_inc, tau_dec, mu;
	int cumu, cumu_rec;
	float rho;
	int p, n, m;
	int M, W;
	real_t *ans;
	float over_relax;
};

void printMatrix(real_t* H, int m, int n);
void read_file(ifstream &Read_other, int * net_vec,
	int * worker_nums, int * worker_ids, int * worker_pos,
	int scale, int M);
void central_cal(string central_file, int M, int W, int m);