#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME AGSTSMO_SENSORLESS
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){
	ssSetNumDiscStates(S, 13); // 8 DISCRETE STATE USED
	if (!ssSetNumInputPorts(S, 1)) return; 
	ssSetInputPortWidth(S, 0, 5); // 6 INPUT
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return; 
	ssSetOutputPortWidth(S, 0, 3); // 10 OUTPUT
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
    Y[0] = X[9]; 
    Y[1] = X[10];
    Y[2] = X[8];
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

    //declare
    real_T k1, k2, k3, k4;
    real_T Va, Vb, Vc, V_alfa, V_beta;
    real_T I_alfa, I_beta;
    real_T surface_alfa, surface_beta, sign_alfa, sign_beta;
    real_T integral_sign_alfa, integral_sign_beta, integral_sign_alfa_old, integral_sign_beta_old;
    real_T I_alfa_estimated_old, I_beta_estimated_old, I_alfa_estimated, I_beta_estimated;
    real_T psi_alfa, psi_beta, phi_alfa, phi_beta;
    real_T integral_phi_alfa, integral_phi_beta, integral_phi_alfa_old, integral_phi_beta_old;
    real_T phi_alfa_old, phi_beta_old;
    real_T integral_psi_alfa, integral_psi_beta, integral_psi_alfa_old, integral_psi_beta_old;
    real_T omega, omega_old, wc, z_alfa, z_beta, psi_alfa_old, psi_beta_old;
    real_T Isd, Isq;
    real_T sqrt23 = 0.816496580927726;
    real_T sqrt32 = 0.866025403784439;
    real_T Tr = Lr / Rr;
    real_T sigma = 1 - (Lm * Lm) / (Ls * Lr);
    real_T gamma = 1 / (sigma * Ls);
    real_T lambda = gamma * Lm / Lr;
    real_T epsilon = Rs / (gamma * Ls);
    wc = 20;
    real_T alpha = wc * dt / (1 + wc * dt);
    real_T beta = 1 / (1 + wc * dt);
    real_T a3 = Lm / (Lr * Ls - Lm * Lm);
    
    // INPUT
    Va = U(0);
    Vb = U(1);
    Vc = U(2);
    I_alfa = U(3);
    I_beta = U(4);

    // STATE
    integral_sign_alfa_old = X[0];
    integral_sign_beta_old = X[1];
    I_alfa_estimated_old = X[2];
    I_beta_estimated_old = X[3];
    integral_psi_alfa_old = X[4];
    integral_psi_beta_old = X[5];
    // integral_phi_alfa_old = X[6];
    // integral_phi_beta_old = X[7];
    phi_alfa_old = X[6];
    phi_beta_old = X[7];
    omega_old = X[8];
    psi_alfa_old = X[11];
    psi_beta_old = X[12];

    // Calculate Alfa-Beta
    V_alfa = sqrt23 * (Va - 0.5 * Vb - 0.5 * Vc);
    V_beta = sqrt23 * (sqrt32 * Vb - sqrt32 * Vc);

    k1 = 200;
    k2 = 1000;
    k3 = 200;
    k4 = 1000;
    // k1 = abs(omega_old) + abs(-1 * (1 / Tr) * phi_alfa_old - omega_old * phi_beta_old + (1 / Tr) * Lm * I_alfa_estimated_old);
    // k2 = 1000;
    // k3 = abs(omega_old) + abs(-1 * (1 / Tr) * phi_beta_old + omega_old * phi_alfa_old + (1 / Tr) * Lm * I_beta_estimated_old);
    // k4 = 1000;

    // CALCULATE SURFACE    
    surface_alfa =  I_alfa_estimated_old - I_alfa;
    surface_beta =  I_beta_estimated_old - I_beta;

    if(surface_alfa > 0){
        sign_alfa = 1;
    }else if(surface_alfa < 0){
        sign_alfa = -1;
    }else{
        sign_alfa = 0;
    }

    if(surface_beta > 0){
        sign_beta = 1; 
    }else if(surface_beta < 0){
        sign_beta = -1;
    }else{
        sign_beta = 0;
    }
    
    // INTEGRAL SIGN   
    integral_sign_alfa = integral_sign_alfa_old + sign_alfa * dt;
    integral_sign_beta = integral_sign_beta_old + sign_beta * dt;

    // CALCULATE PSI    
    psi_alfa = -k1 * pow(fabs(surface_alfa), 0.5) * sign_alfa - k2 * integral_sign_alfa;
    psi_beta = -k3 * pow(fabs(surface_beta), 0.5) * sign_beta - k4 * integral_sign_beta;

    integral_psi_alfa = integral_psi_alfa_old + psi_alfa * dt;
    integral_psi_beta = integral_psi_beta_old + psi_beta * dt;

    if(integral_psi_alfa > 0.5){
        z_alfa = 0.5;
    }else if(integral_psi_alfa < -0.5){
        z_alfa = -0.5;
    }else{
        z_alfa = integral_psi_alfa;
    }

    if(integral_psi_beta > 0.5){
        z_beta = 0.5;
    }else if(integral_psi_beta < -0.5){
        z_beta = -0.5;
    }else{
        z_beta = integral_psi_beta;
    }

    phi_alfa = -integral_psi_alfa;
    phi_beta = -integral_psi_beta;

    I_alfa_estimated = I_alfa_estimated_old + (lambda * psi_alfa - epsilon * I_alfa_estimated_old + gamma * V_alfa) * dt;
    I_beta_estimated = I_beta_estimated_old + (lambda * psi_beta - epsilon * I_beta_estimated_old + gamma * V_beta) * dt;
    omega = omega_old + (0.03 * a3 * (phi_alfa * I_beta_estimated - phi_beta * I_alfa_estimated)) * dt;

    X[0] = integral_sign_alfa;
    X[1] = integral_sign_beta;
    X[2] = I_alfa_estimated;
    X[3] = I_beta_estimated;
    X[4] = integral_psi_alfa;
    X[5] = integral_psi_beta;
    X[6] = phi_alfa;
    X[7] = phi_beta;
    X[8] = omega;
    X[9] = surface_alfa;
    X[10] = surface_beta;
    X[11] = psi_alfa;
    X[12] = psi_beta;
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif
