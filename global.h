using namespace std;

// Fundamental consts
const int n_p = 100; // the number of nodes
const int n_e = n_p-1; // the number of elements
const int d_o_f = 1; //degree of freedom it isn't used here.
const int d_o_s = 6; // degree of SSFEM
const int KL = 3; // order of KL expansion

// Constants for conjugate gradient method
#define Iter_Max 10000 // limit times for iterative methods
#define EPS 1e-14 // error for convergence judgement

// Physical parameters
const double L = 1.0; // length
const double ad = 0.0; // advection velosity
const double Di = 0.01; // diffusion coefficient
const double c_o_v = 0.2; // coefficient of variation
const double stan_len = L;

// Data regarding mesh
vector <double> x; // x-axis coordinate
vector <int> el;   // node number array of elements
vector <int> e_link_start;
vector <int> e_link;
vector <int> n_link_start;
vector <int> n_link; // node link array
vector <int> set_bc_list; // boundary nodes list

// Matrix
vector <double> M; // Mass Matrix
vector <double> M_c; // Correction term Matrix from SUPG method
vector <double> K; // Diffution Matrix 熱伝導行列 力学では剛性行列
vector <double> A; // Advection Matrix
vector <double> A_c; // Correction term Matrix from SUPG method
vector <double> F_body; //熱流束vector　力学では体にかかる外力項
vector <double> F_boundary; //熱流束vector　力学では境界面にかかる外力項
vector <double> u; //solution of this analysis

// Consts for dynamical analysis
const double Courant = 0.1; // When using explicit method, we consider Courant condtion.
const double max_step = 1001;
const double total_time = 1.0;
double dt = total_time/max_step;


//for Eigenvalue problem.
vector <double> Cov;
vector <double> eigenvectors;
vector <double> eigenvalues;

