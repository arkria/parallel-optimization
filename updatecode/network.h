#pragma once
#include<qpOASES.hpp>
#include<string>
#include<fstream>
#include<map>
#include<vector>
#include<set>

USING_NAMESPACE_QPOASES
using namespace std;

class master;
class worker
{
public:
	worker() {};
	worker(int id, int m, int n, float rho_input, int *master_id_input, int *master_pos_input, int worker_tag, ifstream &ReadFile);
	~worker() {};

	// the global rho
	void updatex(float rho);
	void updatelambda(float rho);
	void updateSend(float rho);
	void sendXL(float rho, master* masters);
	float getR_2(master* masters);
	float getR_2();
	float getS_2(float rho);
	float getLambda();
	float getPri();
	float getDual();
	void receiveZ(master *masters, float over_relax);
	
	// for the local rho
	void updatex();
	void updatelambda();
	void updateSend();
	void sendXL(master* masters);
	float getS_2();
	void receiveZ(float over_relax);
	void updateRhP(float R_2, float S_2, float mu, int cumu, float tau_dec, float tau_inc);

	real_t *send_val;
	real_t *z; // update z with over relaxization
	real_t *z_store; //  store the z received from masters
	real_t *x; // primal variable
	real_t *x_temp; // temp of x
	real_t *lambda;
	int m, n; // m is the length of the variable, n is the length of the constraints
	int worker_tag; // local or link worker
	int *master_id; // the id of masters linked with this worker, local worker has only one master, the link one has two
	int *master_pos; //  0 or 1, for the local one, master pos is 0, for the link one, two of its masters are 1 and 0.
	int cumu_rec; // for local rho 
	

private:

	int id;
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
	float rhP;

	/*
	0.5 * x' H x + g x
	s.t. lbA <= Ax <= ubA
	lb < x < ub
	*/
	SQProblem model;

	void updateH(float rho);
	void updateH();
	void updateg(float rho);
	void updateg();
};

class master
{
public:
	master() {};
	master(int m, int worker_num_input, int *con_vec_input, int* worker_pos, ifstream &ReadFile);
	~master() {};
	void updatez(worker *workers);
	void updatez();
	void sendZ(worker *workers); // send information to all the linked workers
	void sendZ(worker *workers, int id); // send information to the specific id workers 
	float getS_2(worker *workers);
	float getS_2();
	float getLambda(worker *workers);
	float getLambda();
	float getPri();
	void printZ();
	real_t **X; // store x received from workers
	real_t **X_temp;
	real_t **L; // store lambda received from workers
	real_t **XL;
	real_t *z; // the consensus variable
	int *worker_pos;

	map<int, int> idToPos;

private:
	int m; // length of the variable
	int *workers_id;
	int worker_num; // number of workers

	real_t *z_temp; // temp of z

};

class ADMM
{
public:
	ADMM(int M, int W, int m, int n, float rho_input,
		int *net_vec_input, int* master_pos, int* worker_nums, int *worker_ids, int* worker_pos,
		float eps_rel, float eps_abs, float tau_inc, float tau_dec, float mu, int cumu, float over_relax,
		string worker_file, string master_file);
	~ADMM() {};
	void solve(int max_step, int output_step, int thread_num);
	void Asyupdate(int max_step, int output_step, int thread_num);
	void AsyFor(int max_step, int output_step, int thread_num);
	real_t * getSolution(real_t *cen_ans);

private:
	master *masters;
	worker *workers;
	float eps_rel, eps_abs, eps_pri, eps_dual, R_2, S_2, Ax, Bz, Alambda;
	float tau_inc, tau_dec, mu;
	int cumu, cumu_rec;
	float rho;
	int p, n, m; // Ax + Bz = c A -> R ^(p x n), m is the length of time domain
	int M, W; // M is the number of masters and W is the number of workers
	real_t *ans;
	float over_relax;

	vector<vector<double>> time_mat;
	vector<double> time_step;
	set<int> arri_set;
	double total_time;

	void updatePri();
	void updatePri2();
	void updatePriFor(int tag, int set_num);
	void updateDual();
	void initial();

	bool check_conv(int tag, int max_step, int output_step);
	bool check_conv(int tag, int max_step, int output_step, float R_2, float S_2, float Ax, float Bz, float Alambda);
	void updateRho();
};

void printMatrix(real_t* H, int m, int n);
void printMatrix(int* H, int m, int n);
void read_file(ifstream &Read_other, int * net_vec, int * master_pos,
	int * worker_nums, int * worker_ids, int * worker_pos,
	int scale, int M);
void central_cal(string central_file, int M, int W, int m, real_t *cen_ans);