#include"network.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include<float.h>
#include<omp.h>
#include<ctime>
USING_NAMESPACE_QPOASES;

using namespace std;
const int local_worker = 1;
const int link_worker = 0;

worker::worker(int m_input, int n_input, float rho_input, 
	int *master_id_input, int worker_tag_input, ifstream &ReadFile)
{
	hot_flag = false;
	m = m_input;
	n = n_input;
	H_temp = new real_t[m*m]();
	g_temp = new real_t[m]();
	A = new real_t[n*m]();
	lbA = new real_t[n]();
	ubA = new real_t[n]();
	lb = new real_t[m]();
	ub = new real_t[m]();
	z = new real_t[m]();
	x = new real_t[m]();
	x_temp = new real_t[m]();

	for (int i = 0; i < m*m; i++)
	{
		ReadFile >> H_temp[i];
	}
	for (int i = 0; i < m; i++)
	{
		ReadFile >> g_temp[i];
	}
	for (int i = 0; i < n*m; i++)
	{
		ReadFile >> A[i];
	}
	for (int i = 0; i < n; i++)
	{
		ReadFile >> ubA[i];
	}
	for (int i = 0; i < n; i++)
	{
		//lbA[i] = FLT_MIN;
		lbA[i] = -10000.0;
	}
	for (int i = 0; i < m; i++)
	{
		ReadFile >> lb[i];
	}
	for (int i = 0; i < m; i++)
	{
		ReadFile >> ub[i];
		//if (i == 0) cout << "worker" << ub[i] << endl;
	}
	for (int i = 0; i < m; i++)
	{
		ReadFile >> x[i];
	}

	SQProblem temp(m, n);
	temp.setPrintLevel(PL_NONE);
	model = temp;
	
	//rho = rho_input;
	lambda = new real_t[m]();
	
	
	send_val = new real_t[m]();
	H = new real_t[m*m]();
	g = new real_t[m]();
	/*updateH();
	updateg();*/
	worker_tag = worker_tag_input;
	master_id = master_id_input;
	/*if (worker_tag == local_worker)
	{
		printf("local worker\n");
		printf("%d\n", master_id[0]);
	}
	else if (worker_tag == link_worker)
	{
		printf("link worker\n");
		printf("%d, %d\n", master_id[0], master_id[1]);
	}*/
}

void worker::updateH(float rho)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) H[i*m + j] = 2 * (H_temp[i*m + j] + 0.5 * rho);
			else H[i*m + j] = 2 * H_temp[i*m + j];
		}
	}
	/*printf("Matrix H is:\n");
	printMatrix(H, m, m);
	printf("Matrix H_temp is:\n");
	printMatrix(H_temp, m, m);*/
}

void worker::updateg(float rho)
{
	for (int i = 0; i < m; i++)
	{
		g[i] = g_temp[i] + lambda[i] - rho * z[i];
	}
}

void worker::updatex(float rho)
{
	updateH(rho), updateg(rho);
	int_t nWSR = m*10;
	if (hot_flag == false)
	{
		hot_flag = true;
		/*printf("\nH is: \n "),	printMatrix(H, m, m);
		printf("\ng is: \n "), printMatrix(g, m, 1);
		printf("\ng_temp is: \n "), printMatrix(g_temp, m, 1);
		printf("\nA is: \n "), printMatrix(A, n, m);
		printf("\nlb is: \n "), printMatrix(lb, m, 1);
		printf("\nub is: \n "), printMatrix(ub, m, 1);
		printf("\nlbA is: \n "), printMatrix(lbA, n, 1);
		printf("\nubA is: \n "), printMatrix(ubA, n, 1);
		printf("\nlambda is: \n "), printMatrix(lambda, m, 1);
		printf("\nz is: \n "), printMatrix(z, m, 1);*/
		model.init(H, g, A, lb, ub, lbA, ubA, nWSR, 0);
		
	}
	else
	{
		model.hotstart(H, g, A, lb, ub, lbA, ubA, nWSR, 0);
	}
	for (int i = 0; i < m; i++)
	{
		x_temp[i] = x[i];
	}
	//x_temp = x;
	model.getPrimalSolution(x);
	/*printf("\n the solved x is: \n "), printMatrix(x, m, 1);
	printf("===========================================================================\n");*/
}

void worker::updatelambda(float rho)
{
	for (int i = 0; i < m; i++)
	{
		lambda[i] = lambda[i] + rho*(x[i] - z[i]);
	}
	/*printf("\n x: \n"), printMatrix(x, m, 1);
	printf("\n z: \n"), printMatrix(z, m, 1);
	printf("\n lambda: \n"), printMatrix(lambda, m, 1);*/
}

//void worker::receivez(master *masters)
//{
//	z = masters[master_id].z;
//}

void worker::updateSend(float rho)
{
	for (int i = 0; i < m; i++)
	{
		send_val[i] = x[i] + lambda[i] / rho;
	}
	/*printf("the send value is:\n");
	printMatrix(x, m, 1);*/
}

float worker::getR_2(master *masters)
{
	float r = 0;
	if (worker_tag == local_worker)
	{
		for (int i = 0; i < m; i++)
		{
			r += pow(x[i] - masters[master_id[0]].z[i], 2);
			//printf("master id is %d, %d %f\n", master_id[0], i, x[i] - masters[master_id[0]].z[i]);
		}
	}
	else if (worker_tag == link_worker)
	{
		for (int i = 0; i < m/2; i++)
		{
			r += pow(x[i] - masters[master_id[0]].z[i], 2);
			//printf("master id is %d, %d %f\n",master_id[0], i, x[i] - masters[master_id[0]].z[i]);
		}
		for (int i = 0; i < m / 2; i++)
		{
			r += pow(x[m/2 + i] - masters[master_id[1]].z[i], 2);
			//printf("master id is %d, %d %f\n", master_id[1], i, x[i] - masters[master_id[0]].z[i]);
		}
	}
	
	return r;
}



float worker::getPri()
{
	float Ax = 0;
	for (int i = 0; i < m; i++)
	{
		Ax = Ax + x[i] * x[i];
	}
	return Ax;
}

float worker::getDual()
{
	float Alambda = 0;
	for (int i = 0; i < m; i++)
	{
		Alambda = Alambda + lambda[i] * lambda[i];
	}
	return Alambda;
}

void worker::receiveZ(master *masters, float over_relax)
{
	if (worker_tag == local_worker)
	{
		
		for (int i = 0; i < m; i++)
		{
			z[i] = over_relax * masters[master_id[0]].z[i] + (1 - over_relax) * x[i];
		}
		/*for (int i = 0; i < m; i++)
		{
			cout << i << " " << z[i] << endl;
		}*/
	}
	else
	{
		//printf("x is:\n");
		//printMatrix(x, m, 1);
		//printf("master id are : %d, %d \n", master_id[0], master_id[1]);
		for (int i = 0; i < m/2; i++)
		{
			z[i] = over_relax * masters[master_id[0]].z[i] + (1 - over_relax) * x[i];
			/*cout << i << " " << z[i] << endl;*/
		}
		//printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		for (int i = 0; i < m/2; i++)
		{
			z[m/2 + i] = over_relax * masters[master_id[1]].z[i] + (1 - over_relax) * x[m/2 + i];
			/*cout << i << " " << z[m+i] << endl;*/
			
			//printf("%d over relax: %f, z: %f, x: %f\n", i+m/2, over_relax, masters[master_id[1]].z[i], x[m/2 + i]);
			
		}
		//printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		/*printf("z is:\n");
		printMatrix(masters[master_id[0]].z, m / 2, 1);
		printMatrix(masters[master_id[1]].z, m / 2, 1);*/
		/*printf("x is:\n");
		printMatrix(x, m, 1);*/
	}
	/*printf("the update value is:\n");
	printMatrix(z, m, 1);*/
}


master::master(int m_input, int worker_num_input, int *con_vec_input, int* worker_pos_input, ifstream &ReadFile)
{
	m = m_input;
	z = new real_t[m]();
	z_temp = z;
	worker_num = worker_num_input;
	workers_id = con_vec_input;
	worker_pos = worker_pos_input;
	/*printf("master\n");
	for (int i = 0; i < worker_num; i++)
	{
		printf("%d ", workers_id[i]);
	}
	printf("\n");*/
}

void master::updatez(worker *workers)
{
	z_temp = z;
	for (int j = 0; j < m; j++)
	{
		z[j] = 0;
	}

	for (int i = 0; i < worker_num; i++)
	{
		//printf("the total worker number is %d\n", worker_num);
		int worker_id = workers_id[i];
		/*cout << "wokrer num " << i << "worker id " << worker_id << endl;*/
		int pos = 0;
		if (worker_pos[i] == 0) pos = 0;
		else if (worker_pos[i] == 1) pos = m;
		for (int j = 0; j < m; j++)
		{
			z[j] = z[j] + workers[worker_id].send_val[pos + j];
			
			//cout <<"z = " << z[j] << " " << workers[worker_id].send_val[j] << endl;
		}
		//printf("worker id: %d, postion: %d\n", worker_id, pos);
	}

	for (int j = 0; j < m; j++)
	{
		z[j] = z[j] / worker_num;
	}
	//printMatrix(z, m, 1);
}

float master::getS_2(worker *workers)
{
	float s = 0;

	for (int i = 0; i < m; i++)
	{
		float temp = 0;
		for (int j = 0; j < worker_num; j++)
		{
			temp += workers[workers_id[j]].x[worker_pos[j] * m + i] - workers[workers_id[j]].x_temp[worker_pos[j] * m + i];
		}
		s += pow(temp, 2);
		//printf("%d %.3f\n", i+1, temp);
	}
	return s;
}

float master::getPri()
{
	float Bz = 0;
	for (int i = 0; i < m; i++)
	{
		Bz = Bz + worker_num * z[i] * z[i];
	}
	return Bz;
}

float master::getLambda(worker *workers)
{
	float s = 0;

	for (int i = 0; i < m; i++)
	{
		float temp = 0;
		for (int j = 0; j < worker_num; j++)
		{
			temp += workers[workers_id[j]].lambda[worker_pos[j] * m + i];
		}
		s += pow(temp, 2);
		//printf("%d %.3f\n", i+1, temp);
	}
	return s;
}

void master::printZ()
{
	printMatrix(z, m, 1);
}

ADMM::ADMM(int M_input, int W_input, int m_input, int n_input, float rho_input, int *net_vec_input, int *worker_nums, int *worker_ids, int* worker_pos, 
	float eps_rel_input, float eps_abs_input, float tau_inc_input, float tau_dec_input, float mu_input, 
	int cumu_input, float over_relax_input, string worker_file, string master_file)
{
	M = M_input, W = W_input;
	p = m_input * M + 2 * m_input * (W - M), n = m_input * M;
	rho = rho_input, eps_rel = eps_rel_input, eps_abs = eps_abs_input;
	R_2 = 0, S_2 = 0;
	tau_inc = tau_inc_input, tau_dec = tau_dec_input, mu = mu_input * mu_input, cumu = cumu_input;
	cumu_rec = 0;
	ans = new real_t[M*m_input]();
	m = m_input;
	over_relax = over_relax_input;
	workers = new worker[W]();
	masters = new master[M]();
	ifstream Read_worker(worker_file);
	ifstream Read_master(master_file);
	for (int i = 0; i < M; i++)
	{
		workers[i] = worker(m_input, n_input, rho, &net_vec_input[i], local_worker, Read_worker);
		if (i == 0) masters[i] = master(m_input, worker_nums[i], &worker_ids[0], &worker_pos[0], Read_master);
		else masters[i] = master(m_input, worker_nums[i] - worker_nums[i-1], &worker_ids[worker_nums[i-1]], &worker_pos[worker_nums[i - 1]], Read_master);
	}
	for (int i = 0; i < (W-M); i++)
	{
		workers[i + M] = worker(2 * m_input, n_input, rho, &net_vec_input[2*i + M], link_worker, Read_worker);
	}
	eps_pri = 0, eps_dual = 0;
}

void ADMM::updatePri()
{
#pragma omp parallel for
	for (int i = 0; i < W; i++)
	{
		//printf("worker %d is being dealed\n", i);
		workers[i].receiveZ(masters, over_relax);
		
		workers[i].updatex(rho);
		workers[i].updatelambda(rho);
		//printf("==========================================================\n");
	}
}

void ADMM::updateDual()
{
#pragma omp parallel for
	for (int i = 0; i < W; i++)
	{
		//printf("update send information, worker %d is being dealed:\n", i);
		workers[i].updateSend(rho);
		/*for (int j = 0; j < workers[i].m; j++)
		{
			printf("worker %d send vector %d is: %f\n", i, j, workers[i].send_val[j]);
		}*/
	}
#pragma omp parallel for
	for (int i = 0; i < M; i++)
	{
		//printf("\nupdate z, master %d is being dealed:\n", i);
		masters[i].updatez(workers);
	}
}

bool ADMM::check_conv(int tag, int max_step, int output_step)
{
	float Ax = 0, Bz = 0, Alambda = 0, pri_val = 0;
	R_2 = 0, S_2 = 0;
#pragma omp parallel for
	for (int i = 0; i < W; i++)
	{
		//printf("worker %d is being dealed with\n", i);
		R_2 += workers[i].getR_2(masters);
		/*cout << i << " " << R_2 << endl;*/
		Ax += workers[i].getPri();
		
	}
#pragma omp parallel for
	for (int i = 0; i < M; i++)
	{
		//printf("master %d is being dealed\n", i + 1);
		S_2 += pow(rho, 2) * masters[i].getS_2(workers);
		Bz += masters[i].getPri();
		Alambda += masters[i].getLambda(workers);
		
	}
	pri_val = Ax > Bz? Ax: Bz;
	eps_pri = (float)(sqrt(p) * eps_abs + sqrt(pri_val) * eps_rel), eps_pri = pow(eps_pri, 2);
	eps_dual = (float)(sqrt(n) * eps_abs + sqrt(Alambda) * eps_rel), eps_dual = pow(eps_dual, 2);
	if (tag % output_step == 0)
	{
		/*cout << "Iteration: " << tag << "  primal criterion: " << sqrt(eps_pri) << "  primal residual: " << sqrt(R_2) 
			<< "  dual criterion: " << sqrt(eps_dual) << "  dual residual: " << sqrt(S_2) << endl;*/
		printf("============= Iteration %d =============\n", tag);
		printf("rho = %f\n", rho);
		printf("eps pri = %f\n", sqrt(eps_pri));
		printf("pri = %f\n", sqrt(R_2));
		printf("eps dual = %f\n", sqrt(eps_dual));
		printf("dual = %f\n", sqrt(S_2));
		/*for (int i = 0; i < M; i++)
		{
			printf("master %d:\n", i + 1);
			masters[i].printZ();
		}*/
		printf("==========================================\n");
	}
	//printf("\n\tstep is %d, pri is: %f, epsi_pri is: %f; dual is %f; epsi_dual is %f\n", tag, sqrt(R_2), sqrt(eps_pri), sqrt(S_2), sqrt(eps_dual));
	if (tag >= max_step || (R_2 < eps_pri && S_2 < eps_dual)) 
		return false;
	else
		return true;

}

void ADMM::updateRho()
{
	if (R_2 > mu * S_2) cumu_rec += 1;
	else if (S_2 > mu * R_2) cumu_rec -= 1;
	
	if (cumu_rec >= cumu)
	{
		rho = rho * tau_inc;
		cumu_rec = 0;
	}
	else if (cumu_rec <= -cumu)
	{
		rho = rho / tau_dec;
		cumu_rec = 0;
	}
}

void ADMM::solve(int max_step, int output_step, int thread_num)
{
	clock_t start, end;
	omp_set_num_threads(thread_num);
	int tag = 0;
	start = clock();
	do {
		updateDual();
		/*for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < m; j++)
			{
				cout << "master " << i << "= " << masters[i].z[j] << endl;
			}
		}*/
		updatePri();
		updateRho();
		tag += 1;
	} while (check_conv(tag, max_step, output_step));
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("The thread number is %d, the total time is %f\n", thread_num, endtime);
}

real_t * ADMM::getSolution()
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < m; j++)
		{
			ans[j + i * m] = masters[i].z[j];
		}
	}
	/*printMatrix(ans, M, m);*/
	return ans;
}

void printMatrix(real_t* H, int m, int n)
{
	for (int j = 0; j < n + 1; j++)
	{
		printf("%10d", j);
	}
	printf("\n");
	for (int i = 0; i < m; i++)
	{
		printf("%10d", i+1);
		for (int j = 0; j < n; j++)
		{
			printf("%10.4f", H[i * n + j]);
		}
		printf("\n");
	}
}

void read_file(ifstream &Read_other, int * net_vec, 
	int * worker_nums, int * worker_ids, int * worker_pos, 
	int scale, int M)
{
	for (int i = 0; i < scale; i++)
	{
		Read_other >> net_vec[i];
	}
	for (int i = 0; i < M; i++)
	{
		Read_other >> worker_nums[i];
	}
	for (int i = 0; i < scale; i++)
	{
		Read_other >> worker_ids[i];
	}
	for (int i = 0; i < scale; i++)
	{
		Read_other >> worker_pos[i];
	}
}

void central_cal(string central_file, int M, int W, int m)
{
	int c_H = M*m;
	int c_A = M * m + (W - M) * 2 * m;
	real_t * H = new real_t[c_H * c_H]();
	real_t * g = new real_t[c_H]();
	real_t * A = new real_t[c_A * c_H]();
	real_t * lbA = new real_t[c_A]();
	real_t * ubA = new real_t[c_A]();
	real_t * lb = new real_t[c_H]();
	real_t * ub = new real_t[c_H]();

	ifstream Read_central(central_file);
	for (int i = 0; i < c_H * c_H; i++)
	{
		Read_central >> H[i];
	}

	for (int i = 0; i < c_H; i++)
	{
		Read_central >> g[i];
	}
	for (int i = 0; i < c_A*c_H; i++)
	{
		Read_central >> A[i];
	}
	for (int i = 0; i < c_A; i++)
	{
		Read_central >> ubA[i];
	}
	for (int i = 0; i < c_A; i++)
	{
		lbA[i] = -10000;
	}
	for (int i = 0; i < c_H; i++)
	{
		Read_central >> lb[i];
	}
	for (int i = 0; i < c_H; i++)
	{
		Read_central >> ub[i];
	}
	/*cout << H[0] << endl;
	cout << g[0] << endl;
	cout << A[0] << endl;
	cout << ubA[0] << endl;
	cout << lb[0] << endl;
	cout << ub[0] << endl;*/
	SQProblem model(c_H, c_A);
	int nWSR = c_H * 10;
	double start, end;
	start = clock();
 	model.setPrintLevel(PL_NONE);
	model.init(H, g, A, lb, ub, lbA, ubA, nWSR, 0);
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("The central time is %f\n", endtime);
}