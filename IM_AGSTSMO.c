#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME IM_AGSTSMO
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
    ssSetNumContStates(S, 5); 
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 4); 
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 4); 
    ssSetNumSampleTimes(S, 1); 

    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); } 

static void mdlInitializeSampleTimes(SimStruct *S) { 
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME); 
    ssSetOffsetTime(S, 0, 0.0); } 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S) { 

    real_T *X0 = ssGetContStates(S); 
    int_T nStates = ssGetNumContStates(S); 
    int_T i; 

    /* initialize the states to 0.0 */ 
    for (i=0; i < nStates; i++) {X0[i] = 0.0;} } 

static void mdlOutputs(SimStruct *S, int_T tid) { 
    real_T *Y = ssGetOutputPortRealSignal(S,0); 
    real_T *X = ssGetContStates(S); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
    real_T ia, ib, ic;
    real_T K    = 0.8164965809; // akar(2/3)
    real_T L    = 0.8660254038; // akar(3/2)
    real_T Ialfas, Ibetas, Ias, Ibs, Ics;
    Ialfas = X[0];
    Ibetas = X[1];

    Ias = K * Ialfas;
    Ibs = K * (-0.5 * Ialfas + L * Ibetas);
    Ics = K * (-0.5 * Ialfas - L * Ibetas);

    Y[0] = Ias;
    Y[1] = Ibs;
    Y[2] = Ics;
    Y[3] = X[4];
} 

#define MDL_DERIVATIVES 
static void mdlDerivatives(SimStruct *S) { 
    real_T *dX = ssGetdX(S); 
    real_T *X = ssGetContStates(S); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
    
    real_T vs_alfa, vs_beta, is_alfa, is_beta, is_alfa_dot, is_beta_dot;
    real_T phir_alfa, phir_beta, phir_alfa_dot, phir_beta_dot, omega, omega_dot;
    real_T TL;
    real_T sqr2_3 = 0.8165;
    real_T sqr3 = 1.732;
    real_T M = 0.315, Ls = 0.337, Lr = 0.337, Rs = 5.4, Rr = 3.8, J = 0.005, p = 2, B = 0.005;

    real_T sigma = 1 - (pow(M,2)/(Ls*Lr));
    real_T Tr = Lr / Rr;
    real_T gamma = 1 / (sigma * Ls);
    real_T lambda = gamma * M / Lr;
    real_T epsilon = Rs / (sigma * Ls);
    real_T m = 1/(sigma*Ls);
    real_T K1 = (1/(sigma*Ls))*(Rs + pow(M,2)/(Tr*Lr));
    real_T K2 = (1/(sigma*Ls))*(M/(Tr*Lr));
    real_T K3 = (1/(sigma*Ls))*(M/Lr);
    real_T K4 = M/Tr;
    real_T K5 = 1/Tr;
    real_T K6 = (M/(J*Lr));


    vs_alfa = sqr2_3*(U(0) - 0.5*U(1) - 0.5*U(2));
    vs_beta = sqr2_3*(0.5*sqr3*U(1) - 0.5*sqr3*U(2));
    TL = U(3);
    
    is_alfa = X[0];
    is_beta = X[1];
    phir_alfa = X[2];
    phir_beta = X[3];
    omega = X[4];
    
    is_alfa_dot = -K1*is_alfa + K2*phir_alfa + K3*p*omega*phir_beta + m*vs_alfa;
    is_beta_dot = -K1*is_beta - K3*p*omega*phir_alfa + K2*phir_beta + m*vs_beta;
    phir_alfa_dot = K4*is_alfa - K5*phir_alfa - p*omega*phir_beta;
    phir_beta_dot = K4*is_beta + p*omega*phir_alfa - K5*phir_beta;
    omega_dot = M/(J*Lr)*(phir_alfa*is_beta - phir_beta*is_alfa) - TL/J-B/J*omega;
    
    dX[0] = is_alfa_dot;
    dX[1] = is_beta_dot;
    dX[2] = phir_alfa_dot;
    dX[3] = phir_beta_dot;
    dX[4] = omega_dot;
} 

static void mdlTerminate(SimStruct *S) 
{} /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /* MEX-file interface mechanism */ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 