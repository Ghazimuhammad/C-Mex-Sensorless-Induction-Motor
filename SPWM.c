#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME SPWM
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){
	if (!ssSetNumInputPorts(S, 1)) return;
	ssSetInputPortWidth(S, 0, 5); // 5 INPUT 
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return;
	ssSetOutputPortWidth(S, 0, 4); // 4 OUTPUT
	ssSetNumSampleTimes(S, 1);
    
	ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); 
} 

static void mdlInitializeSampleTimes(SimStruct *S) { 
	ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME); 
	ssSetOffsetTime(S, 0, 0.0); 
} 

static void mdlOutputs(SimStruct *S, int_T tid) { 
	real_T *Y = ssGetOutputPortRealSignal(S,0); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	real_T t = ssGetT(S);
	
	
    // INPUT
	real_T amp, Fc, Va, Vb, Vc;
    amp = U(0);
    Fc = U(1);
    Va = U(2);
    Vb = U(3);
    Vc = U(4);
    
    real_T halfAmp, modulator;
    int_T i;
    
    // GENERATE TRIANGLE CARRIER
    halfAmp = 0.5 * amp;
    
    modulator = Fc * fmod(t, 1/Fc);
    
    // SPWM ALGORITHM
	for(i = 0; i < 3; i++) {
		if(fabs(U(i + 2)) / halfAmp >= modulator){
			Y[i] = halfAmp; 
		}
		
		else {
			Y[i] = 0.0; 	
		}
		
		if(U(i + 2) < 0.0) { 
			Y[i] = -Y[i]; 
		}
	}
	Y[3] = modulator; 
}

static void mdlTerminate(SimStruct *S) {
	
} /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /* MEX-file interface mechanism */ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif