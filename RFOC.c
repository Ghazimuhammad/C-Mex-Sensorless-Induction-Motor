#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME RFOC
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){
	ssSetNumDiscStates(S, 16); // 8 DISCRETE STATE USED
	if (!ssSetNumInputPorts(S, 1)) return; 
	ssSetInputPortWidth(S, 0, 7); // 6 INPUT
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return; 
	ssSetOutputPortWidth(S, 0, 12); // 10 OUTPUT
	ssSetNumSampleTimes(S, 1); 

	ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 
	| SS_OPTION_DISCRETE_VALUED_OUTPUT));
	
} 

static void mdlInitializeSampleTimes(SimStruct *S){ 
	ssSetSampleTime(S, 0, 1e-4); 
	ssSetOffsetTime(S, 0, 0.0);
} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S){ 
	real_T *X0 = ssGetRealDiscStates(S); 
	int_T nXStates = ssGetNumDiscStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	int_T i; 

	/* initialize the states to 0.0 */ 
	for (i=0; i < nXStates; i++) { 
	X0[i] = 0.0; 
	}
} 

static void mdlOutputs(SimStruct *S, int_T tid){ 
	real_T *Y = ssGetOutputPortRealSignal(S,0); 
	real_T *X = ssGetRealDiscStates(S); 
	
    // OUTPUT
    real_T Va, Vb, Vc;
    real_T Imr, Te, theta_e;
    real_T Isd_ref, Isd_act;
    real_T Isq_ref, Isq_act;
    
    Imr         = X[4];
    theta_e     = X[5];
    Te          = X[6];
    Isd_ref     = X[7];
    Isd_act     = X[8];
    Isq_ref     = X[9];
    Isq_act     = X[10];
    Va          = X[11];
    Vb          = X[12];
    Vc          = X[13];
    
    Y[0] = Va;
    Y[1] = Vb;
    Y[2] = Vc;
    Y[3] = Imr;
    Y[4] = Te;
    Y[5] = theta_e;
    Y[6] = Isd_ref;
    Y[7] = Isd_act;
    Y[8] = Isq_ref;
    Y[9] = Isq_act;
    Y[10] = X[14]; //I alfa
    Y[11] = X[15]; //I beta
} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) {
 
	real_T *X = ssGetRealDiscStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

	real_T dt = 1e-4;
	
	// IM MODEL'S PARAMETER
    real_T N    = 2;
    real_T Lm   = 0.315;
    real_T Ls   = 0.337;
    real_T Lr   = 0.337;
    real_T Rs   = 5.4;
    real_T Rr   = 3.8;
    real_T J    = 0.005;
    real_T B    = 0.0005;
    real_T C    = Lm * Lm - Lr * Ls;
     
    
    // DECLARE VARIABLE
    real_T Isd_ref, Isq_ref;
    
    real_T Ia, Ib, Ic;
    real_T I_alfa, I_beta;
    real_T Isd_act, Isq_act;
    
    real_T Va, Vb, Vc;
    real_T V_alfa, V_beta;
    real_T Vcd, Vcq;
    real_T Vsd, Vsq;
    real_T Vs_max;
    
    real_T error_Isd, error_Isq;
	real_T integral_error_Isd, integral_error_Isq;
	real_T integral_error_Isd_old, integral_error_Isq_old;
    real_T Kidp, Kidi, Kiqp, Kiqi;
    real_T Isd1_ref_new, Isd1_ref_old;
    real_T Isq1_ref_new, Isq1_ref_old;
    real_T Imr_new, Imr_old;
    real_T usd, usq;
    real_T sigma;
    
    real_T omega_m_act;
    
    real_T omega_e;
    real_T theta_e, theta_e_old;
	real_T omega_sl;
    real_T Te;
    
    real_T Tdd  = 0.001;
    real_T Tqd  = 0.001;
    real_T Td   = 0.001;
    real_T T2   = Lr / Rr;
    
    real_T pi = 3.141592654;
    real_T pi2  = 2 * 22 / 7;

    real_T K = 0.8164965809; // akar(2/3)
    real_T L = 0.8660254038; // akar(3/2)
	
    sigma = 1 - Lm * Lm / (Ls * Lr);
    
    Kidp = sigma * Ls / Tdd;
    Kidi = Rs / Tdd;
    Kiqp = sigma * Ls / Tqd;
    Kiqi = Rs / Tqd;
    
    // INPUT
	// Current feedback
	Ia = U(0);
    Ib = U(1);
    Ic = U(2);
    
    // Current reference
    Isd_ref = U(3);
    Isq_ref = U(4);
    
    Vs_max = U(5);
    
    // Speed actual feedback
    omega_m_act = U(6);
    
    // STATE
    integral_error_Isd_old  = X[0];
    integral_error_Isq_old  = X[1];
    Isd1_ref_old            = X[2];
    Isq1_ref_old            = X[3];
    Imr_old                 = X[4];
    theta_e_old             = X[5];
    
    // DECOUPLING    
    Isd1_ref_new    = Isd1_ref_old + ((Isd_ref - Isd1_ref_old) * dt) / Td;
    Isq1_ref_new    = Isq1_ref_old + ((Isq_ref - Isq1_ref_old) * dt) / Td;
    Imr_new         = Imr_old + ((Isd_ref - Imr_old) * dt) / T2;
    
    if(Imr_new <= 0.001) { Imr_new = 0.001; }
        
    // FLUX MODEL
    omega_sl = (Rr * Isq1_ref_new) / (Lr * Isd1_ref_new);
    omega_e = N * omega_m_act + omega_sl;
    theta_e = theta_e_old + omega_e * dt;
    while(theta_e >= pi2) { theta_e -= pi2; }
    while(theta_e < 0) { theta_e += pi2; }
           
    // TRANSFORM ABC TO ALPHA-BETA
	I_alfa = K * (Ia - 0.5 * Ib - 0.5 * Ic);
    I_beta = K * L * (Ib - Ic);
	
	// TRANSFORM ALPHA-BETA TO DQ
	Isd_act = cos(theta_e) * I_alfa + sin(theta_e) * I_beta;
    Isq_act = -sin(theta_e) * I_alfa + cos(theta_e) * I_beta;
    
    Te = N * (1 - sigma) * Ls * Imr_new * Isq_act;
    
    error_Isd = Isd_ref - Isd_act;
    error_Isq = Isq_ref - Isq_act;
    
    integral_error_Isd = integral_error_Isd_old + error_Isd * dt;
    integral_error_Isq = integral_error_Isq_old + error_Isq * dt;
    
    usd = Kidp * error_Isd + Kidi * integral_error_Isd;
    usq = Kiqp * error_Isq + Kiqi * integral_error_Isq;
    
	Vcd = -(omega_e * Ls * sigma * Isq1_ref_new);
    Vcq = omega_e * Ls * sigma * Isd1_ref_new + Ls * (1 - sigma) * omega_e * Imr_new;
    
    Vsd = usd + Vcd;
    Vsq = usq + Vcq;
    
    // TRANSFORM DQ TO ALFA-BETA
    V_alfa = Vsd * cos(theta_e) - Vsq * sin(theta_e);
    V_beta = Vsd * sin(theta_e) + Vsq * cos(theta_e);
    
    // TRANSFORM ALFA-BETA TO ABC  
    Va = K * V_alfa;
    Vb = K * (-0.5 * V_alfa + L * V_beta);
    Vc = K * (-0.5 * V_alfa - L * V_beta);
    
    //Voltage Limiter
    if(Va > Vs_max){
        Va = Vs_max; 
    }
    
    if(Va < (-Vs_max)){ 
        Va = -Vs_max; 
    }
    
    if(Vb > (Vs_max)){ 
        Vb = Vs_max; 
    }
    
    if(Vb < (-Vs_max)){ 
        Vb = -Vs_max; 
    }
    
    if(Vc > Vs_max){ 
        Vc = Vs_max; 
    }
    
    if(Vc < (-Vs_max)){ 
        Vc = -Vs_max; 
    }
    
    
    // UPDATE STATE
    X[0] = integral_error_Isd;
    X[1] = integral_error_Isq;
    X[2] = Isd1_ref_new;
    X[3] = Isq1_ref_new;
    X[4] = Imr_new;
    X[5] = theta_e;
    X[6] = Te;
    X[7] = Isd_ref;
    X[8] = Isd_act;
    X[9] = Isq_ref;
    X[10] = Isq_act;
    X[11] = Va;
    X[12] = Vb;
    X[13] = Vc;
    X[14] = I_alfa;
    X[15] = I_beta;
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif
