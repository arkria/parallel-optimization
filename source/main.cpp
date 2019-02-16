#include<iostream>
#include<qpOASES.hpp>
#include<omp.h>
#include<vector>
#include"network.h"
#include<fstream>
#include<string>

using namespace std;
USING_NAMESPACE_QPOASES;

void main(int argc, char** argv)
{
	string worker_file, master_file, other_file, central_file, fashion;// fashion = 0 (asy) or 1 (syn)
	float rho, tau_inc, tau_dec;
	int ratio, mu, output_step; 
	if (argc == 12)
	{
		worker_file = argv[1], master_file = argv[2], other_file = argv[3], central_file = argv[4], rho = atof(argv[5]);
	    ratio = atoi(argv[6]), tau_inc = atof(argv[7]), tau_dec = atof(argv[8]), mu = atoi(argv[9]), output_step = atoi(argv[10]);
		fashion = argv[11];
	}
	else
	{
		worker_file = "worker.txt", master_file = "master.txt", other_file = "other.txt", central_file = "central.txt", fashion = "adp sy";
		if (fashion == "asy")
		{
			rho = 40, ratio = 50, tau_inc = 1, tau_dec = 1, mu = 10, output_step = 1000;
		}
		else
		{
			rho = 1000, ratio = 1, tau_inc = 2, tau_dec = 2, mu = 11, output_step = 1000;
		}
	} 
	int M = 0, W = 0, m = 15, n = 15, cumu = 22 * ratio; // ratio asy: veh_num = 5, the best ratio is 50
	float eps_rel = 0.0001, eps_abs = 0.0001, over_relax = 1.7;
	int max_step = 100000;//max_step + 1000;
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

	if (fashion == "asy")
		myOpt.AsyFor(max_step, output_step, 100, 1);
	else if (fashion == "adp sy")
		myOpt.synFor(max_step, output_step, 100, W);
	else
		myOpt.solve(max_step, output_step, 100);

	//myOpt.AsyFor(max_step, output_step, 100, 1);
	//myOpt.synFor(max_step, output_step, 100, W);
	//myOpt.Asyupdate(max_step, output_step, 100);
	//myOpt.solve(max_step, output_step, 100);
	//myOpt.Iteration(max_step, output_step, 100, W);
	
	double cen_time = central_cal(central_file, M, W, m, cen_ans);
	myOpt.getSolution(cen_ans, cen_time);
	//system("pause");
}