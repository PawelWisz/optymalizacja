/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include <iostream>
#include <random>
#include"opt_alg.h"
#include"solution.h"




struct Comma final : std::numpunct<char>
{
	char do_decimal_point() const override { return ','; }
};

int main()
{
	try
	{
		cout << "LAB NUMBER " << LAB_NO << endl;
		cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO==0

#elif LAB_NO==1 && LAB_PART==1
		std::uniform_real_distribution<double> unif(-100, 100);
		std::default_random_engine re;
		double alphas[3] = { 2, 5, 10};

		for (int i = 0; i < 3; i++)
		{
			double alpha = alphas[i];
			double output[100][12];
			
			for (int j = 0; j < 100; j++)
			{
				double x0 = unif(re);

				double d = 100, epsilon = 1e-5, gamma = 1e-200; // gamma = d - d_old w lagrange
				int N_max = 1000;

				solution::clear_calls();
				double* p = expansion(x0, d, alpha, N_max);

				solution::clear_calls();
				solution opt_fib = fib(p[0], p[1], epsilon);
				int opt_fib_f_calls = solution::f_calls;

				solution::clear_calls();
				solution opt_lag = lag(p[0], p[1], epsilon, gamma, N_max);
				int opt_lag_f_calls = solution::f_calls;

				output[j][0] = x0;
				output[j][1] = p[0];
				output[j][2] = p[1];
				output[j][3] = p[2];
				output[j][4] = opt_fib.x[0]();
				output[j][5] = opt_fib.y[0]();
				output[j][6] = opt_fib_f_calls;
				output[j][7] = opt_fib.x[0]() >= 62 && opt_fib.x[0]() <= 63 ? 1.0 : 0.0;
				output[j][8] = opt_lag.x[0]();
				output[j][9] = opt_lag.y[0]();
				output[j][10] = opt_lag_f_calls;
				output[j][11] = opt_lag.x[0]() >= 62 && opt_lag.x[0]() <= 63 ? 1.0 : 0.0;;
			}

			std::ofstream out(std::to_string(alpha) + "_lab_1_part_1.csv");

			for (auto& row : output) {
				for (auto col : row)
					out << col << ',';
				out << '\n';
			}
		}
		

		/*
		double output[100][12];
		
		double x0 = 62, d = 70, alpha = 0.5, epsilon = 1e-5, gamma = 1e-200; // gamma = d - d_old w lagrange
		int N_max = 1000;
		double *p = expansion(x0, d, alpha, N_max);

		std::cout << "p[0]: " << p[0] << std::endl;
		std::cout << "p[1]: " <<  p[1] << std::endl;
		std::cout << "p[2]: " <<  p[2] << std::endl;

		solution::clear_calls();
		solution opt_fib = fib(p[0], p[1], epsilon);
		int opt_fib_f_calls = solution::f_calls;
		std::cout << opt_fib << std::endl;

		solution::clear_calls();
		solution opt_lag = lag(p[0], p[1], epsilon, gamma, N_max);
		int opt_lag_f_calls = solution::f_calls;
		std::cout << opt_lag << std::endl;


		output[0][0] = x0;
		output[0][1] = p[0];
		output[0][2] = p[1];
		output[0][3] = p[2];
		output[0][4] = opt_fib.x[0]();
		output[0][5] = opt_fib.y[0]();
		output[0][6] = opt_fib_f_calls;
		output[0][7] = 1.0;
		output[0][8] = opt_lag.x[0]();
		output[0][9] = opt_lag.y[0]();
		output[0][10] = opt_lag_f_calls;
		output[0][11] = 1.0;

		std::ofstream out(std::to_string(alpha) + "_lab_1_part_1.csv");

		for (auto& row : output) {
			for (auto col : row)
				out << col << ',';
			out << '\n';
		}
		*/
		
#elif LAB_NO==1 && LAB_PART==2
		double epsilon = 1e-5, gamma = 1e-200;
		int N_max = 1000;

		matrix ab_F(1, 1, 200);
		solution::clear_calls();
		solution opt_f = fib(-100, 100, epsilon, &ab_F);
		std::cout << opt_f << std::endl;
		std::cout << solution::f_calls << std::endl;
		std::cout << ab_F << std::endl;
		
		matrix ab_L(1, 1, 200);
		solution::clear_calls();
		solution opt_l = lag(-100, 100, epsilon, gamma, N_max, &ab_L);
		std::cout << opt_l << std::endl;
		std::cout << ab_L << std::endl;

#elif LAB_NO==1 && LAB_PART==3
		double epsilon = 1e-5, gamma = 1e-200;
		int N_max = 1000;
		
		int row = 1001;
		int col = 3;

		matrix ab_F(row, col);
		solution::clear_calls();
		solution opt_f = fib(0.0001, 0.01, epsilon, &ab_F);
		std::cout << opt_f << std::endl;
		
		std::ofstream out_f("lab_1_part_3_fib.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_f << ab_F(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_f << ab_F(i, j) << "\n";
				}
		}

		matrix ab_L(row, col);
		solution::clear_calls();
		solution opt_l = lag(0.0001, 0.01, epsilon, gamma, N_max, &ab_L);
		std::cout << opt_l << std::endl;

		std::ofstream out_l("lab_1_part_3_lag.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_l << ab_L(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_l << ab_L(i, j) << "\n";
				}
		}
		
#elif LAB_NO==2 && LAB_PART==1
		double epsilon = 1e-3, alpha_HJ = 0.5, aplha_R = 2, beta = 0.5;
		int Nmax = 5000;

		double steps[3] = { 0.1, 0.5, 1.0 };
		auto isGlobalMin = [epsilon](double x1, double x2) {
			x1 = abs(x1);
			x2 = abs(x2);
			return (x1 >= 0 && x1 < epsilon) && (x2 >= 0 && x2 < epsilon);
		};

		for (int i = 0; i < 3; i++)
		{
			double s = steps[i];
			matrix s0 = matrix(2, 1, s);
			double output[100][12];

			for (int j = 0; j < 100; j++)
			{
				matrix x0 = 2 * rand_mat(2, 1) - 1;

				solution::clear_calls();
				solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax);
				int opt_hj_f_calls = solution::f_calls;

				solution::clear_calls();
				solution opt_R = Rosen(x0, s0, aplha_R, beta, epsilon, Nmax);
				int opt_rosen_f_calls = solution::f_calls;

				output[j][0] = x0(0);
				output[j][1] = x0(1);
				output[j][2] = opt_HJ.x(0);
				output[j][3] = opt_HJ.x(1);
				output[j][4] = opt_HJ.y(0);
				output[j][5] = opt_hj_f_calls;
				output[j][6] = isGlobalMin(opt_HJ.x(0), opt_HJ.x(1)) ? 1.0 : 0;
				output[j][7] = opt_R.x(0);
				output[j][8] = opt_R.x(1);
				output[j][9] = opt_R.y(0);
				output[j][10] = opt_rosen_f_calls;
				output[j][11] = isGlobalMin(opt_R.x(0), opt_R.x(1)) ? 1.0 : 0;
			}

			std::ofstream out(std::to_string(s) + "_lab_2_part_1.csv");

			for (auto& row : output) {
				for (auto col : row)
					out << col << ',';
				out << '\n';
			}
		}
		
#elif LAB_NO==2 && LAB_PART==2
		double s = 0.1, epsilon = 1e-3, alpha_HJ = 0.5, aplha_R = 2, beta = 0.5;
		int Nmax = 5000;
		matrix s0 = matrix(2, 1, s);
		matrix x0 = matrix(2, 1);
		x0(0) = -0.0432063;
		x0(1) = -0.343887;
		cout << x0 << endl << endl;
		
		solution::clear_calls();
		matrix ud_HJ(1, 2);
		solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax, &ud_HJ);
		cout << opt_HJ << endl << endl;
		cout << ud_HJ << endl << endl;
		
		solution::clear_calls();
		matrix ud_R(1, 2);
		solution opt_R = Rosen(x0, s0, aplha_R, beta, epsilon, Nmax, &ud_R);
		cout << opt_R << endl;
		cout << ud_R << endl;

#elif LAB_NO==2 && LAB_PART==3
		double s = 0.1, epsilon = 1e-3, alpha_HJ = 0.5, aplha_R = 2, beta = 0.5;
		int Nmax = 5000;
		int row = 1001;
		int col = 2;
		matrix s0 = matrix(2, 1, s);
		matrix x0 = matrix(2, new double[2]{ 1, 1 });
		cout << x0 << endl << endl;

		solution::clear_calls();
		matrix ud_HJ(row, col);
		solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax, &ud_HJ);
		cout << opt_HJ << endl << endl;

		std::ofstream out_hj("lab_2_part_3_hj.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_hj << ud_HJ(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_hj << ud_HJ(i, j) << "\n";
				}
		}

		solution::clear_calls();
		matrix ud_R(row, col);
		solution opt_R = Rosen(x0, s0, aplha_R, beta, epsilon, Nmax, &ud_R);
		cout << opt_R << endl;

		std::ofstream out_r("lab_2_part_3_r.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_r << ud_R(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_r << ud_R(i, j) << "\n";
				}
		}
#elif LAB_NO==3 && LAB_PART==1

std::ofstream pen_zew;
pen_zew.imbue(std::locale(std::locale::classic(), new Comma));
std::ofstream pen_wew;
pen_wew.imbue(std::locale(std::locale::classic(), new Comma));
pen_zew.open("Lab03_P01_F_zew.csv");
pen_wew.open("Lab03_P01_F_wew.csv");

 for (int i = 0; i < 100; i++) {
	matrix x0, a = 4; // a znajduje siê w g3


	// zrobiæ losowanie tak, ¿eby by³o w obaszarze dopuszczalnym
	do {
		x0 = 4 * rand_mat(2, 1) + 1;
	} while (norm(x0) > a(0));

	pen_zew << x0(0) << ";" << x0(1);

	//cout << x0 << endl << endl;


	// zewnetrzna
	double c0 = 1, dc = 2, epsilon = 1e-5;

	int Nmax = 10000;

	// 4,49

	solution opt_zew = pen(x0, c0, dc, epsilon, Nmax, &a);
	double result = sqrt(pow(opt_zew.x(0), 2) + pow(opt_zew.x(1), 2));
	cout << opt_zew << endl << endl;
	//cout << sqrt(pow(opt_zew.x(0), 2) + pow(opt_zew.x(1), 2))  << endl << endl; // 4.000001
	cout << result << endl;

	pen_zew << ";" << opt_zew.x(0) << ";" << opt_zew.x(1) << ";" << result << ";" << opt_zew.y(0) << ";" << opt_zew.f_calls << "\n";


	solution::clear_calls();


	// wewnetrzna
	c0 = 10;
	dc = 0.5;

	solution opt_wew = pen(x0, c0, dc, epsilon, Nmax, &a);
	cout << opt_wew << endl << endl;
	//cout << sqrt(pow(opt_wew.x(0), 2) + pow(opt_zew.x(1), 2)) << endl << endl; // 4.000001
	result = sqrt(pow(opt_wew.x(0), 2) + pow(opt_wew.x(1), 2));
	cout << result << endl;

	pen_wew << ";" << opt_wew.x(0) << ";" << opt_wew.x(1) << ";" << result << ";" << opt_wew.y(0) << ";" << opt_wew.f_calls << "\n";

	solution::clear_calls();
}

pen_zew.close();
pen_wew.close();

#elif LAB_NO==3 && LAB_PART==2

matrix x0(2, 1, 2), c = 1;
solution test(x0);
test.fit_fun(nullptr, &c);
cout << test << endl;

#elif LAB_NO==4 && LAB_PART==1

std::ofstream SD_method;
std::ofstream CG_method;
std::ofstream Newton_method;
SD_method.imbue(std::locale(std::locale::classic(), new Comma));
CG_method.imbue(std::locale(std::locale::classic(), new Comma));
Newton_method.imbue(std::locale(std::locale::classic(), new Comma));
SD_method.open("Lab04_P01_SD.csv");
CG_method.open("Lab04_P01_CG.csv");
Newton_method.open("Lab04_P01_Newton.csv");

for (int i = 0; i < 100; i++) {
	matrix x0 = 20 * rand_mat(2, 1) - 10;
	double epsilon = 1e-3;
	//double h0 = 0.05;  //h0 ustala krok, jesli minus to jest zmiennokrokowe 
	double h0 = 0.12;
	//double h0 = -0.05;
	int Nmax = 10000;

	solution opt;
	opt = SD(x0, h0, epsilon, Nmax);

	cout << "Metoda najszybszego spadku" << endl;
	cout << opt << endl << endl;
	SD_method << x0(0) << ";" << x0(1) << ";" << opt.x(0) << ";" << opt.x(1) << ";" << opt.y(0) << ";" << opt.f_calls << ";" << opt.g_calls << "\n";
	solution::clear_calls();

	opt = CG(x0, h0, epsilon, Nmax);

	cout << "Metoda gradientow sprzezonych" << endl;
	cout << opt << endl << endl;
	CG_method << opt.x(0) << ";" << opt.x(1) << ";" << opt.y(0) << ";" << opt.f_calls << ";" << opt.g_calls << "\n";
	solution::clear_calls();

	opt = Newton(x0, h0, epsilon, Nmax);  // if h0 < 0 g_calls i h_calls musi byc pomiedzy 2-3

	cout << "Metoda Newtona" << endl;
	cout << opt << endl << endl;
	Newton_method << opt.x(0) << ";" << opt.x(1) << ";" << opt.y(0) << ";" << opt.f_calls << ";" << opt.g_calls << ";" << opt.H_calls << "\n";
	solution::clear_calls();
}
SD_method.close();
CG_method.close();
Newton_method.close();

#elif LAB_NO==4 && LAB_PART==2

std::ofstream SD_method;
std::ofstream CG_method;
std::ofstream Newton_method;
SD_method.imbue(std::locale(std::locale::classic(), new Comma));
CG_method.imbue(std::locale(std::locale::classic(), new Comma));
Newton_method.imbue(std::locale(std::locale::classic(), new Comma));
SD_method.open("Lab04_P02_SD.csv");
CG_method.open("Lab04_P02_CG.csv");
Newton_method.open("Lab04_P02_Newton.csv");


matrix x0(2, new double[2]{ 1.61514 , -1.16231});
double epsilon = 1e-3;
//double h0 = 0.05;
//double h0 = 0.12;
double h0 = -0.05;
int row=100;
int col=2;
int Nmax = 10000;

matrix ud_SD(trans(x0));
solution opt;
cout << "Metoda najszybszego spadku" << endl << endl;
opt = SD(x0, h0, epsilon, Nmax, &ud_SD);
cout << opt << endl;
cout << "maciez ud" << endl;
cout << ud_SD << endl;
SD_method << ud_SD << endl;
solution::clear_calls();

matrix ud_CG(trans(x0));
cout << "Metoda gradientow sprzezonych" << endl << endl;
opt = CG(x0, h0, epsilon, Nmax, &ud_CG);
cout << opt << endl;
cout << "maciez ud" << endl;
cout << ud_CG << endl;
CG_method << ud_CG << endl;
solution::clear_calls();


matrix ud_Newton(trans(x0));
cout << "Metoda Newtona" << endl << endl;
opt = Newton(x0, h0, epsilon, Nmax, &ud_Newton);
cout << opt << endl;
cout << "maciez ud" << endl;
cout << ud_Newton << endl;
Newton_method << ud_Newton << endl;
solution::clear_calls();

SD_method.close();
CG_method.close();
Newton_method.close();

#elif LAB_NO==4 && LAB_PART==3

double steps[3] = { 0.01, 0.001, 0.0001 };
double epsilon = 1e-5;
int Nmax = 10000;

for (int i = 0; i < 3; i++)
{
	matrix x0(3, new double[3]{ 0, 0, 0 });
	double h0 = steps[i];

	solution opt = CG(x0, h0, epsilon, Nmax);
	cout << opt << endl;
	solution::clear_calls();
}

#elif LAB_NO==5 && LAB_PART==1

#elif LAB_NO==5 && LAB_PART==2

#endif
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
