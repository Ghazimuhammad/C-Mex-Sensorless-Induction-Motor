#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME IP
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
	ssSetNumDiscStates(S, 3); // 3 DISCRETE STATE USED
	if (!ssSetNumInputPorts(S, 1)) return; 
	ssSetInputPortWidth(S, 0, 4); // 5 INPUT
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return; 
	ssSetOutputPortWidth(S, 0, 2); // 2 OUTPUT
	ssSetNumSampleTimes(S, 1); 

	ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 
	| SS_OPTION_DISCRETE_VALUED_OUTPUT));
} 

static void mdlInitializeSampleTimes(SimStruct *S){ 
	ssSetSampleTime(S, 0, 1e-3); 
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

    real_T Tref, error;
    
    Tref = X[1];
    error = X[2];
    
	Y[0] = Tref;
    
    Y[1] = error;
} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) {
	real_T *X = ssGetRealDiscStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

	real_T dt = 1e-3;

    // INPUT
    real_T Kp, Ki, speed_ref, Isd_ref, speed_act;
    Kp          = U(0);
	Ki          = U(1);
    speed_ref   = U(2);
	speed_act   = U(3);
    
    
    real_T Tmax = 1.4;
	real_T error;
	real_T integral, integral_old;
    real_T Tref;


	// STATE
	integral_old = X[0];
	
    // IP ALGORITHM
	error       = speed_ref - speed_act; 
    integral    = integral_old + error * dt; 
    Tref        = Ki * integral - Kp * speed_act;
    
    if(Tref >= Tmax) {
        Tref = Tmax;
        integral = integral_old;
    }
    
    else if(Tref <= -Tmax) {
        Tref = -Tmax;
        integral = integral_old;
    }
    
	X[0] = integral;
	X[1] = Tref;	
    X[2] = error;
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif
