// Number of  bins
const int N_Q2  = 3;
const int N_Nu  = 3;
const int N_Zh = 4;
const int N_Pt2 = 6;
const int N_Phi = 6;
const float Delta_PT2 = 3.0/N_Pt2;
const float Delta_Phi = 360.0/N_Phi;
const int N_PION = 3;
// Number of solids targets
const int N_STARGETS = 3;

// Limits
const float Zh_MIN = 0.0;
const float Zh_MAX = 1.0;
const float Pt2_MIN = 0.0;
const float Pt2_MAX = 3.0;


// Binning
const float Q2_BINS[N_Q2+1] = {1, 1.3, 1.8, 4.0};
const float Nu_BINS[N_Nu+1] = {2.2, 3.2, 3.7, 4.26};
const float Zh_BINS[N_Zh+1] = {0, 0.2, 0.4, 0.6, 1.};
const float Pt2_BINS[N_Pt2+1] = {0.0, 0.5, 1, 1.5, 2 ,2.5, 3};
const float Phi_BINS[N_Phi+1] = {-180, -120, -60, 0, 60, 120, 180};
