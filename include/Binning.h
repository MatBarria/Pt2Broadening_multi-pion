// Number of  bins
const int N_Q2  = 3;
const int N_Nu  = 3;
const int N_Zh = 8;
const int N_Pt2 = 80;
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
const float Zh_BINS[N_Zh+1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.};
//const float Pt2_BINS[N_Pt2+1] = {0.0, 0.5, 1, 1.5, 2 ,2.5, 3};
const float Phi_BINS[N_Phi+1] = {-180, -120, -60, 0, 60, 120, 180};
const float Pt2_BINS[N_Pt2+1] = {0, 0.0375, 0.075, 0.1125, 0.15, 0.1875, 0.225, 0.2625, 0.3, 0.3375, 0.375,
                                0.4125, 0.45, 0.4875, 0.525, 0.5625, 0.6, 0.6375, 0.675, 0.7125, 0.75, 0.7875,
                                0.825, 0.8625, 0.9, 0.9375, 0.975, 1.0125, 1.05, 1.0875, 1.125, 1.1625, 1.2,
                                1.2375, 1.275, 1.3125, 1.35, 1.3875, 1.425, 1.4625, 1.5, 1.5375, 1.575, 1.6125,
                                1.65, 1.6875, 1.725, 1.7625, 1.8, 1.8375, 1.875, 1.9125, 1.95, 1.9875, 2.025,
                                2.0625, 2.1, 2.1375, 2.175, 2.2125, 2.25, 2.2875, 2.325, 2.3625, 2.4, 2.4375,
                                2.475, 2.5125, 2.55, 2.5875, 2.625, 2.6625, 2.7, 2.7375, 2.775, 2.8125, 2.85,
                                2.8875, 2.925, 2.9625, 3};
