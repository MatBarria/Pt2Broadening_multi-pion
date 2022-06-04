// Number of  bins
const int N_Q2  = 3;
const int N_Nu  = 3;
const int N_Zh = 4;
const int N_Pt2 = 90;
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
const int ZH_SUM = 1;


// Binning
const float Q2_BINS[N_Q2+1] = {1, 1.3, 1.8, 4.0};
const float Nu_BINS[N_Nu+1] = {2.2, 3.2, 3.7, 4.26};
const float Zh_BINS[N_Zh+1] = {0, 0.2, 0.4, 0.6 ,1.0};
//const float Pt2_BINS[N_Pt2+1] = {0.0, 0.5, 1, 1.5, 2 ,2.5, 3};
const float Phi_BINS[N_Phi+1] = {-180, -120, -60, 0, 60, 120, 180};
const float Pt2_BINS[N_Pt2+1] = {0, 0.0333333, 0.0666667, 0.1, 0.133333, 0.166667, 0.2, 0.233333, 0.266667, 0.3, 0.333333,
                                0.366667, 0.4, 0.433333, 0.466667, 0.5, 0.533333, 0.566667, 0.6, 0.633333, 0.666667, 0.7,
                                0.733333, 0.766667, 0.8, 0.833333, 0.866667, 0.9, 0.933333, 0.966667, 1, 1.03333, 1.06667,
                                1.1, 1.13333, 1.16667, 1.2, 1.23333, 1.26667, 1.3, 1.33333, 1.36667, 1.4, 1.43333,
                                1.46667, 1.5, 1.53333, 1.56667, 1.6, 1.63333, 1.66667, 1.7, 1.73333, 1.76667, 1.8,
                                1.83333, 1.86667, 1.9, 1.93333, 1.96667, 2, 2.03333, 2.06667, 2.1, 2.13333, 2.16667,
                                2.2, 2.23333, 2.26667, 2.3, 2.33333, 2.36667, 2.4, 2.43333, 2.46667, 2.5, 2.53333,
                                2.56667, 2.6, 2.63333, 2.66667, 2.7, 2.73333, 2.76667, 2.8, 2.83333, 2.86667, 2.9,
                                2.93333, 2.96667, 3};
