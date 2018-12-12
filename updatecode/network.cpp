#include"network.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include<float.h>
#include<omp.h>
#include<ctime>
#include<map>
#include<algorithm>
USING_NAMESPACE_QPOASES;

using namespace std;
/*========================= other function and data structure ======================================================*/
void printMatrix(real_t* H, int m, int n)
{
	for (int j = 0; j < n + 1; j++)
	{
		printf("%10d", j);
	}
	printf("\n");
	for (int i = 0; i < m; i++)
	{
		printf("%10d", i + 1);
		for (int j = 0; j < n; j++)
		{
			printf("%10.4f", H[i * n + j]);
		}
		printf("\n");
	}
}

void printMatrix(int* H, int m, int n)
{
	for (int j = 0; j < n + 1; j++)
	{
		printf("%10d", j);
	}
	printf("\n");
	for (int i = 0; i < m; i++)
	{
		printf("%10d", i + 1);
		for (int j = 0; j < n; j++)
		{
			printf("%10d", H[i * n + j]);
		}
		printf("\n");
	}
}

void read_file(ifstream &Read_other, int * net_vec, int * master_pos,
	int * worker_nums, int * worker_ids, int * worker_pos,
	int scale, int M)
{
	for (int i = 0; i < scale; i++)
	{
		Read_other >> net_vec[i];
	}

	for (int i = 0; i < scale; i++)
	{
		Read_other >> master_pos[i];
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

void central_cal(string central_file, int M, int W, int m, real_t *cen_ans)
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
	model.getPrimalSolution(cen_ans);
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("The central time is %f\n", endtime);
}

struct compare_temp
{
	int id;
	double time;
};

bool mycompare(compare_temp c1, compare_temp c2)
{
	return(c1.time < c2.time);
}
const int local_worker = 1;
const int link_worker = 0;

/*======================== worker ===================================================*/
worker::worker(int id_input, int m_input, int n_input, float rho_input,
	int *master_id_input, int *master_pos_input, int worker_tag_input, ifstream &ReadFile)
{
	id = id_input;
	hot_flag = false;
	m = m_input;
	n = n_input;
	rhP = rho_input;
	cumu_rec = 0;
	H_temp = new real_t[m*m]();
	g_temp = new real_t[m]();
	A = new real_t[n*m]();
	lbA = new real_t[n]();
	ubA = new real_t[n]();
	lb = new real_t[m]();
	ub = new real_t[m]();
	z = new real_t[m]();
	z_store = new real_t[m]();
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
		lbA[i] = -10000.0;
	}
	for (int i = 0; i < m; i++)
	{
		ReadFile >> lb[i];
	}
	for (int i = 0; i < m; i++)
	{
		ReadFile >> ub[i];
	}
	for (int i = 0; i < m; i++)
	{
		ReadFile >> x[i];
	}

	SQProblem temp(m, n);
	temp.setPrintLevel(PL_NONE);
	model = temp;
	lambda = new real_t[m]();

	send_val = new real_t[m]();
	H = new real_t[m*m]();
	g = new real_t[m]();
	worker_tag = worker_tag_input;
	master_id = master_id_input;
	master_pos = master_pos_input;
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
}

void worker::updateH()
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) H[i*m + j] = 2 * (H_temp[i*m + j] + 0.5 * rhP);
			else H[i*m + j] = 2 * H_temp[i*m + j];
		}
	}
}

void worker::updateg(float rho)
{
	for (int i = 0; i < m; i++)
	{
		g[i] = g_temp[i] + lambda[i] - rho * z[i];
	}
}

void worker::updateg()
{
	for (int i = 0; i < m; i++)
	{
		g[i] = g_temp[i] + lambda[i] - rhP * z[i];
	}
}

void worker::updatex(float rho)
{
	updateH(rho), updateg(rho);
	int_t nWSR = m * 10;
	if (hot_flag == false)
	{
		hot_flag = true;
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
	model.getPrimalSolution(x);
}

void worker::updatex()
{
	updateH(), updateg();
	int_t nWSR = m * 10;
	if (hot_flag == false)
	{
		hot_flag = true;
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
	model.getPrimalSolution(x);
}

void worker::updatelambda(float rho)
{
	for (int i = 0; i < m; i++)
	{
		lambda[i] = lambda[i] + rho*(x[i] - z[i]);
	}
}

void worker::updatelambda()
{
	for (int i = 0; i < m; i++)
	{
		lambda[i] = lambda[i] + rhP*(x[i] - z[i]);
	}
}

void worker::updateSend(float rho)
{
	for (int i = 0; i < m; i++)
	{
		send_val[i] = x[i] + lambda[i] / rho;
	}
	/*printf("the send value is:\n");
	printMatrix(x, m, 1);*/
}

void worker::updateSend()
{
	for (int i = 0; i < m; i++)
	{
		send_val[i] = x[i] + lambda[i] / rhP;
	}
	/*printf("the send value is:\n");
	printMatrix(x, m, 1);*/
}

void worker::sendXL(float rho, master* masters)
{
	for (int i = 0; i < m; i++)
	{
		send_val[i] = x[i] + lambda[i] / rho;
	}

	if (worker_tag == local_worker)
	{
		int pos = masters[master_id[0]].idToPos[id];
		//printf("pos =: %d and worker id is: %d\n", pos, id);
		for (int i = 0; i < m; i++)
		{
			masters[master_id[0]].X[pos][i] = x[i];
			masters[master_id[0]].X_temp[pos][i] = x_temp[i];
			masters[master_id[0]].L[pos][i] = lambda[i];
			masters[master_id[0]].XL[pos][i] = send_val[i];
		}
	}
	else
	{
		int pos = masters[master_id[0]].idToPos[id];
		//printf("pos =: %d, and workerid is : %d, start pos is: %d and %d\n", pos, id, master_pos[0] * m / 2, master_pos[1] * m / 2);
		for (int i = 0; i < m / 2; i++)
		{
			masters[master_id[0]].X[pos][i] = x[master_pos[0] * m / 2 + i];
			masters[master_id[0]].X_temp[pos][i] = x_temp[master_pos[0] * m / 2 + i];
			masters[master_id[0]].L[pos][i] = lambda[master_pos[0] * m / 2 + i];
			masters[master_id[0]].XL[pos][i] = send_val[master_pos[0] * m / 2 + i];
		}
		pos = masters[master_id[1]].idToPos[id];
		for (int i = 0; i < m / 2; i++)
		{
			masters[master_id[1]].X[pos][i] = x[master_pos[1] * m / 2 + i];
			masters[master_id[1]].X_temp[pos][i] = x_temp[master_pos[1] * m / 2 + i];
			masters[master_id[1]].L[pos][i] = lambda[master_pos[1] * m / 2 + i];
			masters[master_id[1]].XL[pos][i] = send_val[master_pos[1] * m / 2 + i];
		}
	}

}

void worker::sendXL(master* masters)
{
	for (int i = 0; i < m; i++)
	{
		send_val[i] = x[i] + lambda[i] / rhP;
	}

	if (worker_tag == local_worker)
	{
		int pos = masters[master_id[0]].idToPos[id];
		//printf("pos =: %d and worker id is: %d\n", pos, id);
		for (int i = 0; i < m; i++)
		{
			masters[master_id[0]].X[pos][i] = x[i];
			masters[master_id[0]].X_temp[pos][i] = x_temp[i];
			masters[master_id[0]].L[pos][i] = lambda[i];
			masters[master_id[0]].XL[pos][i] = send_val[i];
		}
	}
	else
	{
		int pos = masters[master_id[0]].idToPos[id];
		//printf("pos =: %d, and workerid is : %d, start pos is: %d and %d\n", pos, id, master_pos[0] * m / 2, master_pos[1] * m / 2);
		for (int i = 0; i < m / 2; i++)
		{
			masters[master_id[0]].X[pos][i] = x[master_pos[0] * m / 2 + i];
			masters[master_id[0]].X_temp[pos][i] = x_temp[master_pos[0] * m / 2 + i];
			masters[master_id[0]].L[pos][i] = lambda[master_pos[0] * m / 2 + i];
			masters[master_id[0]].XL[pos][i] = send_val[master_pos[0] * m / 2 + i];
		}
		pos = masters[master_id[1]].idToPos[id];
		for (int i = 0; i < m / 2; i++)
		{
			masters[master_id[1]].X[pos][i] = x[master_pos[1] * m / 2 + i];
			masters[master_id[1]].X_temp[pos][i] = x_temp[master_pos[1] * m / 2 + i];
			masters[master_id[1]].L[pos][i] = lambda[master_pos[1] * m / 2 + i];
			masters[master_id[1]].XL[pos][i] = send_val[master_pos[1] * m / 2 + i];
		}
	}

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
		for (int i = 0; i < m / 2; i++)
		{
			r += pow(x[i] - masters[master_id[0]].z[i], 2);
			//printf("master id is %d, %d %f\n",master_id[0], i, x[i] - masters[master_id[0]].z[i]);
		}
		for (int i = 0; i < m / 2; i++)
		{
			r += pow(x[m / 2 + i] - masters[master_id[1]].z[i], 2);
			//printf("master id is %d, %d %f\n", master_id[1], i, x[i] - masters[master_id[0]].z[i]);
		}
	}

	return r;
}

float worker::getR_2()
{
	float r = 0;
	for (int i = 0; i < m; i++)
	{
		r += pow(x[i] - z_store[i], 2);
	}
	return r;
}

float worker::getS_2()
{
	float s = 0;

	for (int i = 0; i < m; i++)
	{
		s += pow(rhP * (x[i] - x_temp[i]), 2);
	}
	return s;
}

float worker::getS_2(float rho)
{
	float s = 0;

	for (int i = 0; i < m; i++)
	{
		s += pow(rho * (x[i] - x_temp[i]), 2);
	}
	return s;
}

float worker::getLambda()
{
	float s = 0;

	for (int i = 0; i < m; i++)
	{
		s += pow(lambda[i], 2);
	}
	return s;
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
	}
	else
	{
		for (int i = 0; i < m / 2; i++)
		{
			z[i] = over_relax * masters[master_id[0]].z[i] + (1 - over_relax) * x[i];
			/*cout << i << " " << z[i] << endl;*/
		}
		for (int i = 0; i < m / 2; i++)
		{
			z[m / 2 + i] = over_relax * masters[master_id[1]].z[i] + (1 - over_relax) * x[m / 2 + i];
		}
	}
}

void worker::receiveZ(float over_relax)
{
	for (int i = 0; i < m; i++)
	{
		z[i] = over_relax * z_store[i] + (1 - over_relax) * x[i];
	}
}

void worker::updateRhP(float R_2, float S_2, float mu, int cumu, float tau_dec, float tau_inc)
{
	if (R_2 > mu * S_2) cumu_rec += 1;
	else if (S_2 > mu * R_2) cumu_rec -= 1;

	if (cumu_rec >= cumu)
	{
		rhP = rhP * tau_inc;
		cumu_rec = 0;
	}
	else if (cumu_rec <= -cumu)
	{
		rhP = rhP / tau_dec;
		cumu_rec = 0;
	}
}

/*========================= master =====================================================*/
master::master(int m_input, int worker_num_input, int *con_vec_input, int* worker_pos_input, ifstream &ReadFile)
{
	m = m_input;
	z = new real_t[m]();
	z_temp = z;
	worker_num = worker_num_input;
	workers_id = con_vec_input;
	worker_pos = worker_pos_input;
	for (int i = 0; i < worker_num; i++)
	{
		idToPos[workers_id[i]] = i;
	}
	X = new real_t*[worker_num]();
	for (int i = 0; i < worker_num; i++)
	{
		X[i] = new real_t[m]();
	}

	X_temp = new real_t*[worker_num]();
	for (int i = 0; i < worker_num; i++)
	{
		X_temp[i] = new real_t[m]();
	}

	L = new real_t*[worker_num]();
	for (int i = 0; i < worker_num; i++)
	{
		L[i] = new real_t[m]();
	}

	XL = new real_t*[worker_num]();
	for (int i = 0; i < worker_num; i++)
	{
		XL[i] = new real_t[m]();
	}
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
		/*printf("worker id is: %d\n", i);
		printf("XL:\n");
		printMatrix(XL[idToPos[workers_id[i]]], 1, m);
		printf("primal:\n");
		printMatrix(&workers[worker_id].send_val[pos], 1, m);*/
		for (int j = 0; j < m; j++)
		{
			z[j] = z[j] + workers[worker_id].send_val[pos + j];
		}
	}

	for (int j = 0; j < m; j++)
	{
		z[j] = z[j] / worker_num;
	}
	//printMatrix(z, m, 1);
}

void master::updatez()
{
	z_temp = z;
	for (int j = 0; j < m; j++)
	{
		z[j] = 0;
	}

	for (int i = 0; i < worker_num; i++)
	{
		for (int j = 0; j < m; j++)
		{
			z[j] = z[j] + XL[idToPos[workers_id[i]]][j];
		}
	}

	for (int j = 0; j < m; j++)
	{
		z[j] = z[j] / worker_num;
	}
}

void master::sendZ(worker *workers)
{
	//printMatrix(workers_id, 1, worker_num);
	for (int i = 0; i < worker_num; i++)
	{
		//printf("worker id is %d: \n", workers_id[i]);
		for (int j = 0; j < m; j++)
		{
			if (workers[workers_id[i]].worker_tag == local_worker)
			{
				workers[workers_id[i]].z_store[j] = z[j];
			}
			else
			{
				workers[workers_id[i]].z_store[worker_pos[i] * m + j] = z[j];
			}

		}
	}
}

void master::sendZ(worker *workers, int id)
{
	//printMatrix(workers_id, 1, worker_num);
	for (int i = 0; i < worker_num; i++)
	{
		if (id == workers_id[i])
		{
			for (int j = 0; j < m; j++)
			{
				if (workers[workers_id[i]].worker_tag == local_worker)
				{
					workers[workers_id[i]].z_store[j] = z[j];
				}
				else
				{
					workers[workers_id[i]].z_store[worker_pos[i] * m + j] = z[j];
				}

			}
		}
		//printf("worker id is %d: \n", workers_id[i]);

	}
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

float master::getS_2()
{
	float s = 0;

	for (int i = 0; i < m; i++)
	{
		float temp = 0;
		for (int j = 0; j < worker_num; j++)
		{
			temp += X[j][i] - X_temp[j][i];
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

float master::getLambda()
{
	float s = 0;

	for (int i = 0; i < m; i++)
	{
		float temp = 0;
		for (int j = 0; j < worker_num; j++)
		{
			temp += L[j][i];
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

/*======================== ADMM ===================================================================================*/

ADMM::ADMM(int M_input, int W_input, int m_input, int n_input, float rho_input, int *net_vec_input, int *master_pos, int *worker_nums, int *worker_ids, int* worker_pos,
	float eps_rel_input, float eps_abs_input, float tau_inc_input, float tau_dec_input, float mu_input,
	int cumu_input, float over_relax_input, string worker_file, string master_file)
{
	M = M_input, W = W_input;
	p = m_input * M + 2 * m_input * (W - M), n = m_input * M;
	rho = rho_input, eps_rel = eps_rel_input, eps_abs = eps_abs_input;
	R_2 = 0, S_2 = 0, Ax = 0, Bz = 0, Alambda = 0;
	tau_inc = tau_inc_input, tau_dec = tau_dec_input, mu = mu_input * mu_input, cumu = cumu_input;
	cumu_rec = 0;
	ans = new real_t[M*m_input]();
	m = m_input;
	over_relax = over_relax_input;
	workers = new worker[W]();
	masters = new master[M]();
	total_time = 0;

	for (int i = 0; i < W; i++)
	{
		arri_set.insert(i);
	}

	ifstream Read_worker(worker_file);
	ifstream Read_master(master_file);
	for (int i = 0; i < M; i++)
	{
		workers[i] = worker(i, m_input, n_input, rho, &net_vec_input[i], &master_pos[i], local_worker, Read_worker);
		if (i == 0) masters[i] = master(m_input, worker_nums[i], &worker_ids[0], &worker_pos[0], Read_master);
		else masters[i] = master(m_input, worker_nums[i] - worker_nums[i - 1], &worker_ids[worker_nums[i - 1]], &worker_pos[worker_nums[i - 1]], Read_master);
		/*printMatrix(masters[i].workers_id, 1, masters[i].worker_num);*/
	}
	for (int i = 0; i < (W - M); i++)
	{
		workers[i + M] = worker(i + M, 2 * m_input, n_input, rho, &net_vec_input[2 * i + M], &master_pos[2 * i + M], link_worker, Read_worker);
	}
	eps_pri = 0, eps_dual = 0;
}

void ADMM::updatePri()
{
	R_2 = 0, Ax = 0, S_2 = 0, Alambda = 0;
#pragma omp parallel for
	for (int i = 0; i < W; i++)
	{
		workers[i].receiveZ(over_relax);
		workers[i].updatex(rho);
		workers[i].updatelambda(rho);
		R_2 += workers[i].getR_2();
		/*cout << i << " " << R_2 << endl;*/
		Ax += workers[i].getPri();
		S_2 += workers[i].getS_2(rho);
		Alambda += workers[i].getLambda();
		workers[i].sendXL(rho, masters);
		//printf("==========================================================\n");
	}
}

void ADMM::updatePri2()
{
	float r_2 = 0, s_2 = 0;
	R_2 = 0, Ax = 0, S_2 = 0, Alambda = 0;
#pragma omp parallel for
	for (int i = 0; i < W; i++)
	{
		workers[i].receiveZ(over_relax);
		workers[i].updatex();
		workers[i].updatelambda();
		r_2 = workers[i].getR_2();
		R_2 += r_2;
		/*cout << i << " " << R_2 << endl;*/
		Ax += workers[i].getPri();
		s_2 = workers[i].getS_2();
		S_2 += s_2;
		Alambda += workers[i].getLambda();
		workers[i].sendXL(masters);
		workers[i].updateRhP(r_2, s_2, mu, cumu, tau_dec, tau_inc);
		//printf("==========================================================\n");
	}
}

void ADMM::updatePriFor(int tag, int set_num)
{
	R_2 = 0, Ax = 0, S_2 = 0, Alambda = 0;
	clock_t t1, t2;
	vector<double> temp(W, 0);
	for (int i = 0; i < W; i++)
	{
		if (arri_set.find(i) != arri_set.end())
		{
			t1 = clock();
			workers[i].receiveZ(over_relax);
			workers[i].updatex(rho);
			workers[i].updatelambda(rho);
			workers[i].sendXL(rho, masters);
			t2 = clock();
			temp[i] = (double)(t2 - t1) / CLOCKS_PER_SEC;
		}
		else
		{
			temp[i] = time_mat[tag - 1][i] - time_step[tag - 1];
		}
		//printf("==========================================================\n");
	}
	for (int i = 0; i < W; i++)
	{
		if (arri_set.find(i) != arri_set.end())
		{
			R_2 += workers[i].getR_2();
			/*cout << i << " " << R_2 << endl;*/
			Ax += workers[i].getPri();
			S_2 += workers[i].getS_2(rho);
			Alambda += workers[i].getLambda();
		}
		//printf("==========================================================\n");
	}

	time_mat.push_back(temp);
	vector<compare_temp> sort_time;
	for (int i = 0; i < W; i++)
	{
		sort_time.push_back({ i, temp[i] });
	}
	sort(sort_time.begin(), sort_time.end(), mycompare);
	arri_set.clear();
	for (int i = 0; i < set_num; i++)
	{
		arri_set.insert(sort_time[i].id);
	}
	int i = set_num - 1;
	while (i < W - 1 && abs(sort_time[i].time - sort_time[i + 1].time) < 1e-6)
	{
		arri_set.insert(sort_time[i + 1].id);
		i += 1;
	}
	time_step.push_back(sort_time[set_num - 1].time);
	total_time += sort_time[set_num - 1].time;
}

void ADMM::updateDual()
{
	Bz = 0;
#pragma omp parallel for
	for (int i = 0; i < M; i++)
	{
		masters[i].updatez();
		Bz += masters[i].getPri();
		masters[i].sendZ(workers);
	}
}

void ADMM::initial()
{
#pragma omp parallel for
	for (int i = 0; i < W; i++)
	{
		workers[i].sendXL(rho, masters);
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
	pri_val = Ax > Bz ? Ax : Bz;
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
		printf("The arrival set number is %d and the element are:\n", arri_set.size());
		for (set<int>::iterator it = arri_set.begin(); it != arri_set.end(); it++)
		{
			printf("%d ", *it);
		}
		printf("\n");
		printf("==========================================\n");
	}
	//printf("\n\tstep is %d, pri is: %f, epsi_pri is: %f; dual is %f; epsi_dual is %f\n", tag, sqrt(R_2), sqrt(eps_pri), sqrt(S_2), sqrt(eps_dual));
	if (tag >= max_step || (R_2 < eps_pri && S_2 < eps_dual))
		return false;
	else
		return true;

}

bool ADMM::check_conv(int tag, int max_step, int output_step, float R_2, float S_2, float Ax, float Bz, float Alambda)
{
	float pri_val = 0;
	pri_val = Ax > Bz ? Ax : Bz;
	eps_pri = (float)(sqrt(p) * eps_abs + sqrt(pri_val) * eps_rel), eps_pri = pow(eps_pri, 2);
	eps_dual = (float)(sqrt(n) * eps_abs + sqrt(Alambda) * eps_rel), eps_dual = pow(eps_dual, 2);
	//eps_dual = (float)(sqrt(Alambda) * eps_rel), eps_dual = pow(eps_dual, 2);
	if (tag % output_step == 0)
	{
		printf("============= Iteration %d =============\n", tag);
		printf("rho = %f\n", rho);
		printf("eps pri = %f\n", sqrt(eps_pri));
		printf("pri = %f\n", sqrt(R_2));
		printf("eps dual = %f\n", sqrt(eps_dual));
		printf("dual = %f\n", sqrt(S_2));
		printf("==========================================\n");
	}
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

real_t * ADMM::getSolution(real_t *cen_ans)
{
	float max_error = 0.0;
	float first_error = 0.0;
	float temp = 0;
	int pos_x = 0;
	int pos_y = 0;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < m; j++)
		{
			ans[j + i * m] = masters[i].z[j];
			temp = abs(ans[j + i * m] - cen_ans[j + i * m]) / 0.5236;
			/*max_error = max_error > temp ? max_error : temp;*/
			if (j == 0)
			{
				first_error = first_error > temp ? first_error : temp;
			}
			if (max_error < temp)
			{
				max_error = temp;
				pos_x = i;
				pos_y = j;
			}
		}
	}
	printf("The max error of the first step is %f%%\n", first_error * 100);
	printf("The max error is %f%%, the position is (%d, %d)\n", max_error * 100, pos_x + 1, pos_y + 1);

	/*printMatrix(cen_ans, M, m);
	printMatrix(ans, M, m);*/
	return ans;
}

// ******************** the solution fashion you can edit ***************************
void ADMM::solve(int max_step, int output_step, int thread_num)
{
	clock_t start, end;
	omp_set_num_threads(thread_num);
	int tag = 0;
	start = clock();
	initial();
	do {
		updateDual();
		updatePri();
		updateRho();
		tag += 1;
	} while (check_conv(tag, max_step, output_step));
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("The total time is %f, the iteration number is %d\n", endtime, tag);
}

void ADMM::Asyupdate(int max_step, int output_step, int thread_num)
{
	clock_t start, end;
	omp_set_num_threads(thread_num);
	int tag = 0;
	start = clock();
	int flag = 1;

	initial();
	updateDual();

#pragma omp parallel for
	for (int i = 0; i < W; i++)
	{
		//printf("worker %d is being dealed\n", i);
		do {
			tag++;
			//printf("the step is %d, and the thread is %d\n", tag, i);
			workers[i].receiveZ(over_relax);
			workers[i].updatex(rho);
			workers[i].updatelambda(rho);
			workers[i].sendXL(rho, masters);

#pragma omp critical
			{
				if (workers[i].worker_tag == local_worker)
				{
					masters[workers[i].master_id[0]].updatez();
					masters[workers[i].master_id[0]].sendZ(workers, i);
				}
				else
				{
					masters[workers[i].master_id[0]].updatez();
					masters[workers[i].master_id[0]].sendZ(workers, i);
					masters[workers[i].master_id[1]].updatez();
					masters[workers[i].master_id[1]].sendZ(workers, i);
				}
			}

		} while (tag <= max_step);

		//printf("==========================================================\n");

	}
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("The total time is %f, the iteration number is %d\n", endtime, tag);
}

void ADMM::AsyFor(int max_step, int output_step, int thread_num)
{
	clock_t start, end;
	omp_set_num_threads(thread_num);
	int tag = 0;
	start = clock();
	initial();
	do {
		updateDual();
		updatePriFor(tag, 1);
		updateRho();
		tag += 1;
	} while (check_conv(tag, max_step, output_step));
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("The total time is %f, the iteration number is %d\n", total_time, tag);
}