#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME IM
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
	ssSetNumContStates(S, 6); // 6 CONTINUOUS STATE
	if (!ssSetNumInputPorts(S, 1)) return; 
	ssSetInputPortWidth(S, 0, 4); // 4 INPUT
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return; 
	ssSetOutputPortWidth(S, 0, 8); // 8 OUTPUT
	ssSetNumSampleTimes(S, 1); 
	ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); 
} 

static void mdlInitializeSampleTimes(SimStruct *S) {
	ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME); 
	ssSetOffsetTime(S, 0, 0.0);
} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S) { 

	real_T *X0 = ssGetContStates(S); 
	int_T nStates = ssGetNumContStates(S); 
	int_T i; 

	/* initialize the states to 0.0 */ 
	for (i=0; i < nStates; i++) {
		X0[i] = 0.0;
	}
} 

static void mdlOutputs(SimStruct *S, int_T tid) { 
	real_T *Y = ssGetOutputPortRealSignal(S,0); 
	real_T *X = ssGetContStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	
	// IM MODEL'S PARAMETER
    real_T N    = 2;
    real_T Lm   = 0.2279;
    real_T Ls   = 0.2349;
    real_T Lr   = 0.2349;
    real_T Rs   = 2.76;
    real_T Rr   = 2.9;
    real_T J    = 0.0436;
    real_T B    = 0.0005;
    real_T C    = Lm * Lm - Lr * Ls;
    
    real_T Ialfas_dot, Ibetas_dot;
    real_T Ialfar_dot, Ibetar_dot;
    real_T Ialfas, Ibetas;
    real_T Ialfar, Ibetar;
    real_T omega_m, omega_m_dot;
    real_T theta_m, theta_m_dot;
    real_T theta_e, theta_e_dot;  
 
    real_T psi_alfa, psi_beta;
    real_T Te, TL;
    real_T Imr;
    
    real_T Vas, Vbs, Vcs;
    real_T Valfas, Vbetas;
    real_T Ias, Ibs, Ics;
    
    real_T K    = 0.8164965809; // akar(2/3)
    real_T L    = 0.8660254038; // akar(3/2)
    real_T pi2  = 2 * 22 / 7;
    
    
    // STATE VARIABLE
    Ialfas = X[0];
    Ibetas = X[1];
    Ialfar = X[2];
    Ibetar = X[3];
    omega_m = X[4];
    theta_m = X[5];
    
    
    // TRANSFORM ALPHA-BETA TO ABC
    Ias = K * Ialfas;
    Ibs = K * (-0.5 * Ialfas + L * Ibetas);
    Ics = K * (-0.5 * Ialfas - L * Ibetas);
    
    
    // OUTPUT
    psi_alfa = Lr * Ialfar + Lm * Ialfas;
    psi_beta = Lr * Ibetar + Lm * Ibetas;
    
    Imr = (sqrt(psi_alfa * psi_alfa + psi_beta * psi_beta)) / Lm;
    
    Te = N * Lm * (Ibetas * Ialfar - Ialfas * Ibetar);
    
    if(psi_beta == 0) {
        psi_beta = 1e-6;
    }
    
    if(psi_alfa < 0.0 && psi_beta < 0.0) {
        theta_e = pi2 + atan2(psi_beta, psi_alfa);
    }
    
    else if(psi_alfa < 0.0 && psi_beta > 0.0) {
        theta_e = atan2(psi_beta, psi_alfa);
    }
    
    else if(psi_alfa >= 0 && psi_beta < 0.0) {
        theta_e = pi2 + atan2(psi_beta, psi_alfa);
    }
    
    else if(psi_alfa >= 0 && psi_beta > 0.0) {
        theta_e = atan2(psi_beta, psi_alfa);
    }
    
    while(theta_m > pi2) { theta_m -= pi2; }
    while(theta_m < 0.0) { theta_m += pi2; }
    
    // Current phase A
    Y[0] = Ias;
    
    // Current phase B
    Y[1] = Ibs;
    
    // Current phase C
    Y[2] = Ics;
    
    // Magnetizing Rotor Current
    Y[3] = Imr;
    
    // Electric Torque
    Y[4] = Te;
    
    // Rotor flux angle
    Y[5] = theta_e;
    
    // Rotor angle
    Y[6] = theta_m;

    // Rotor speed
    Y[7] = omega_m;
} 

#define MDL_DERIVATIVES 
static void mdlDerivatives(SimStruct *S) { 
	
	real_T *dX = ssGetdX(S); 
	real_T *X = ssGetContStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	
	// IM MODEL'S PARAMETER
    real_T N    = 2;
    real_T Lm   = 0.2279;
    real_T Ls   = 0.2349;
    real_T Lr   = 0.2349;
    real_T Rs   = 2.76;
    real_T Rr   = 2.9;
    real_T J    = 0.0436;
    real_T B    = 0.0005;
    real_T C    = Lm * Lm - Lr * Ls;
    
    real_T Ialfas_dot, Ibetas_dot;
    real_T Ialfar_dot, Ibetar_dot;
    real_T Ialfas, Ibetas;
    real_T Ialfar, Ibetar;
    real_T omega_m, omega_m_dot;
    real_T theta_m, theta_m_dot;
    real_T theta_e, theta_e_dot;  
 
    real_T psi_alfa, psi_beta;
    real_T Te, TL;
    real_T Imr;
    
    real_T Vas, Vbs, Vcs;
    real_T Valfas, Vbetas;
    real_T Ias, Ibs, Ics;
    
    real_T K    = 0.8164965809; // akar(2/3)
    real_T L    = 0.8660254038; // akar(3/2)
    real_T pi2  = 2 * 22 / 7;
    
    
    // INPUT
	Vas = U(0);
    Vbs = U(1);
    Vcs = U(2);
    TL = U(3);
    
  
    // STATE VARIABLE
    Ialfas = X[0];
    Ibetas = X[1];
    Ialfar = X[2];
    Ibetar = X[3];
    omega_m = X[4];
    theta_m = X[5];
    
    
    // TRANSFORM
    //ABC TO ALPHA BETA
    Valfas = K * (Vas - 0.5 * Vbs - 0.5 * Vcs);
    Vbetas = K * L * (Vbs - Vcs);
    
    
    // IM MODEL
    // State Equation
    Ialfas_dot = (Rs * Lr * Ialfas - N * omega_m * Lm * Lm * Ibetas - Rr * Lm * Ialfar - N * omega_m * Lr * Lm * Ibetar - Lr * Valfas) / C;
    Ibetas_dot = (N * omega_m * Lm * Lm * Ialfas + Rs * Lr * Ibetas + N * omega_m * Lr * Lm * Ialfar - Rr * Lm * Ibetar - Lr * Vbetas) / C;
    Ialfar_dot = -(Rs * Lm * Ialfas - N * omega_m * Lm * Ls * Ibetas - Rr * Ls * Ialfar - N * omega_m * Lr * Ls * Ibetar - Lm * Valfas) / C;
    Ibetar_dot = -(N * omega_m * Lm * Ls * Ialfas + Rs * Lm * Ibetas + N * omega_m * Lr * Ls * Ialfar - Rr * Ls * Ibetar - Lm * Vbetas) / C;
    
    Te = N * Lm * (Ibetas * Ialfar - Ialfas * Ibetar);
    omega_m_dot = (Te - TL - B * omega_m) / J;
    theta_m_dot = omega_m;
    

    // STATE DERIVATIVES
    dX[0] = Ialfas_dot;
    dX[1] = Ibetas_dot;
    dX[2] = Ialfar_dot;
    dX[3] = Ibetar_dot;
    dX[4] = omega_m_dot;
    dX[5] = theta_m_dot;
} 

static void mdlTerminate(SimStruct *S) 
{} /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /* MEX-file interface mechanism */ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 