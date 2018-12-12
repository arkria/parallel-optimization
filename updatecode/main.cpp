#include<iostream>
#include<qpOASES.hpp>
#include<omp.h>
#include<vector>
#include"network.h"
#include<fstream>

using namespace std;
USING_NAMESPACE_QPOASES;

void main()
{
	string worker_file = "worker.txt", master_file = "master.txt", other_file = "other.txt";
	int M = 0, W = 0, m = 15, n = 15, cumu = 22;
	float rho = 1000.0, eps_rel = 0.0001, eps_abs = 0.0001, tau_inc = 2, tau_dec = 2, mu = 10, over_relax = 1.7;
	int max_step = 50000, output_step = 1000;
	ifstream Read_other(other_file);
	Read_other >> M, Read_other >> W;
	int scale = M + 2 * (W - M);
	int * net_vec = new int[scale]();
	int * master_pos = new int[scale]();
	int * worker_nums = new int[M]();
	int * worker_ids = new int[scale]();
	int * worker_pos = new int[scale]();

	real_t * cen_ans = new real_t[M*m]();

	read_file(Read_other, net_vec, master_pos, worker_nums, worker_ids, worker_pos, scale, M);

	/*int net_vec[14] = { 0, 1, 2, 3, 0, 1, 0, 3, 1, 2, 1, 3, 2, 3 };
	int worker_nums[4] = { 3, 7, 10, 14 };
	int worker_ids[14] = {0, 4, 5, 1, 4, 6, 7, 2, 6, 8, 3, 5, 7, 8};
	int worker_pos[14] = { 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1 };*/ 
	/*real_t H[2 * 3] = { 0,1,2,3,4,5 };
	printMatrix(H, 2, 3);*/

	ADMM myOpt(M, W, m, n, rho, net_vec, master_pos, worker_nums, worker_ids, worker_pos, eps_rel,
		eps_abs, tau_inc, tau_dec, mu, cumu, over_relax, worker_file, master_file);
	printf("The vehicle number is %d\n", M);
	//myOpt.AsyFor(max_step, output_step, 100);
	//myOpt.Asyupdate(max_step, output_step, 100);
	myOpt.solve(max_step, output_step, 100);
	central_cal("central.txt", M, W, m, cen_ans);
	myOpt.getSolution(cen_ans);
	system("pause");
}