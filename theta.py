# theta.py -- A catalysis reaction modeled in Python.
# Copyright (c) Philip Conrad 2013 -- MIT License


from mpmath import *


#----------------------------------------------------------
# REACTION RATE CONSTANTS:

k_f = [ '6.249e7',  '1.831e3',  '9.567e2',  '8.447e3', '1.863e5', '5.509e8',   '5.982', '2.106e10',  '7.139e4',  '2.243e8',
        '2.418e7',  '1.247e8',  '1.100e2', '5.791e12', '1.114e9', '9.955e3', '5.396e2',  '3.706e3',  '2.705e8',  '7.040e9',
        '5.501e8',  '2.335e4', '1.630e10',  '6.622e2', '6.464e2', '2.109e8',   '8.910',  '3.268e5',  '1.890e5', '9.014e11',
        '7.631e3',  '6.303e2',  '1.075e2',  '9.362e7', '9.540e4', '2.636e8', '3.368e8', '1.615e10', '3.290e-3',  '1.004e3',
        '1.457e5',  '2.380e2',  '3.845e7',  '3.778e7', '9.432e3', '1.666e3', '3.094e8',  '1.557e7',  '6.575e1',  '1.372e2',
          '3.003', '3.044e13', '1.047e16', '2.092e12', '1.020e3',
        '5.056e9',  '1.396e9',  '4.831e9',  '9.712e6',   '4.000', '2.403e6', '1.404e-1' #interesting values left on their own row
]

k_b = [ '5.885e4', '3.070e6',  '2.885e7', '8.560e1',    '8.721',  '3.131e5', '6.828e-12',  '4.823e3',   '1.020e4',    '1.566',
        '5.024e5', '3.056e6', '6.318e-1', '1.247e8', '8.518e-3', '6.388e10',   '1.625e7',  '1.514e5',   '4.864e5',  '1.941e2',
        '1.750e6', '8.974e1', '5.505e-2', '1.555e1',  '2.380e7',    '1.986',   '1.865e7',    '1.668',   '1.108e7',    '1.962',
       '1.902e11', '5.235e1',  '1.311e7',   '2.729',  '2.606e8',  '5.388e6',     '4.689',  '1.170e4', '2.457e-15', '1.833e-6',
       '7.225e-5', '5.142e8', '4.046e12', '9.921e9',  '3.620e7',  '2.431e7',  '1.802e13',  '2.232e8',   '7.117e7',  '6.635e7',
        '7.879e7', '1.230e8',  '1.740e8', '1.640e8',  '6.696e7',
        '2.126e8',   '1.611',  '1.848e6', '1.828e2',  '3.558e8', '1.593e5', '1.336e10'  #interesting values left on their own row
]


#----------------------------------------------------------
# SYSTEM VARIABLE VALUES:

initialState = [
     '3.9e-10',  '4.2e-11', '3.92e-14',  '2.61e-3',  '4.75e-6', '6.21e-12',  '5.01e-8', '2.13e-11',  '5.32e-9',  '5.66e-2',
     '2.80e-7',  '2.60e-4',  '4.99e-8',  '6.69e-2',  '4.50e-3',    '0.554', '3.02e-17', '1.08e-18',  '2.21e-3',  '2.59e-6',
     '1.22e-6', '6.07e-19', '5.81e-15', '9.66e-16', '1.72e-08', '7.49e-20', '2.92e-10', '1.29e-11', '4.37e-11', '1.22e-15', 
    '1.85e-10', '4.04e-15', '6.26e-13', '2.53e-11', '3.77e-17', '3.60e-09', '7.43e-10', '1.25e-12', '4.31e-06'
]

minBounds = [
    '1e-10', '1e-30', '1e-30', '1e-40', '1e-20', '1e-25', '1e-10', '1e-20', '1e-20', '1e-30',
    '1e-20', '1e-20', '1e-20', '1e-20', '1e-20', '1e-30', '1e-25', '1e-20', '1e-20', '1e-20',
    '1e-20', '1e-25', '1e-20', '1e-20', '1e-20', '1e-20', '1e-20', '1e-20', '1e-20', '1e-30',
    '1e-20', '1e-20', '1e-20', '1e-20', '1e-30', '1e-20', '1e-20', '1e-20', '1e-30'
]

maxBounds = [
    '0.80', '0.01', '0.01', '0.80', '0.01', '0.01', '0.50', '0.01', '0.01', '0.80',
    '0.01', '0.90', '0.01',  '0.8', '0.80', '0.80', '0.80', '0.80', '0.80', '0.50',
    '0.01', '0.01', '0.01', '0.10', '0.80', '0.50', '0.10', '0.10', '0.50', '0.10',
    '0.10', '0.10', '0.10', '0.10', '0.10', '0.10', '0.10', '0.10', '0.80'
]

#Partial pressures, Boltzman constant, and temperature:
p_A     = mpf('1.0')          #guaiacol
p_H2    = mpf('1.0')          #H2
p_M     = mpf('1.0e-7')       #catechol
p_S     = mpf('1.0e-7')       #phenol
p_X     = mpf('1.0e-7')       #benzene
p_K     = mpf('1.0e-7')       #anisole
p_CO    = mpf('1.0e-7')
p_H2O   = mpf('0.0')
p_CH4   = mpf('0.0')
p_CH3OH = mpf('0.0')
KB      = mpf('8.6173324e-5') #Boltzman constant
T       = mpf('573')          #temperature


#----------------------------------------------------------
# FUNCTIONS:

#initList :: [pyFloat] -> [mpf]
def initList(src):
    return list(map(lambda x: mpf(x), src))


#----------------------------------------------------------
# MAIN:

if __name__ == '__main__':
    #--------------------------------------
    # SETUP:

    #Set precision to be very high:
    mp.dps = 50

    #Pretty print our results:
    mp.pretty = True

    #Initialize all the lists of constants:
    k_f             = initList(k_f)
    k_b             = initList(k_b)
    initialState    = initList(initialState)
    minBounds       = initList(minBounds)
    maxBounds       = initList(maxBounds)

    #Globals allowing us to examine the reaction rates and steady state approximation later:
    fr = []
    rr = []
    tr = []
    ss = []


    #--------------------------------------
    # SOLVER:

    def coverage(theta_A,
                 theta_E,
                 theta_F,
                 theta_G,
                 theta_I,
                 theta_L,
                 theta_M,
                 theta_N,
                 theta_O,
                 theta_R,
                 theta_S,
                 theta_U,
                 theta_V,
                 theta_X,
                 free,
                 theta_H,
                 theta_OH,
                 theta_H2O,
                 theta_CH,
                 theta_CH2,
                 theta_CH3,
                 theta_CH4,
                 theta_CHO,
                 theta_CH2O,
                 theta_CH3O,
                 theta_CH3OH,
                 theta_W,
                 theta_T,
                 theta_B,
                 theta_J,
                 theta_P,
                 theta_C,
                 theta_D,
                 theta_K,
                 theta_Q,
                 theta_Y,
                 theta_Z,
                 theta_Z1,
                 theta_CO
                ):

        forwardRates = [
            lambda : k_f[0]*p_A*free*free*free*free,
            lambda : k_f[1]*theta_A*theta_H,
            lambda : k_f[2]*theta_A*theta_H,        #special
            lambda : k_f[3] *theta_A*free,          #special
            lambda : k_f[4]*theta_A*free,
            lambda : k_f[5]*theta_A*free,
            lambda : k_f[6]*theta_A*free,
            lambda : k_f[7]*theta_A*free,
            lambda : k_f[8]*theta_B*free*free,
            lambda : k_f[9] *theta_C*free,          #special
            lambda : k_f[10]*theta_D*theta_H,       #special
            lambda : k_f[11]*theta_E*theta_H,
            lambda : k_f[12]*theta_F*free,
            lambda : k_f[13]*theta_F*free,
            lambda : k_f[14]*theta_F,
            lambda : k_f[15]*theta_G*theta_H,
            lambda : k_f[16]*theta_G*free,
            lambda : k_f[17]*theta_I*free,
            lambda : k_f[18]*theta_I*free*free,
            lambda : k_f[19]*theta_J,
            lambda : k_f[20]*theta_K*free,          #special
            lambda : k_f[21]*theta_L*free,
            lambda : k_f[22]*theta_L*free,
            lambda : k_f[23]*theta_M*free,
            lambda : k_f[24]*theta_M*theta_H,
            lambda : k_f[25]*theta_N*theta_H,
            lambda : k_f[26]*theta_O*free,
            lambda : k_f[27]*theta_O*free,
            lambda : k_f[28]*theta_P*theta_H,
            lambda : k_f[29]*theta_Q*free,          #special
            lambda : k_f[30]*theta_R*theta_H,
            lambda : k_f[31]*theta_S*free,
            lambda : k_f[32]*theta_S*theta_H,
            lambda : k_f[33]*theta_T*free,
            lambda : k_f[34]*theta_U*theta_H,
            lambda : k_f[35]*theta_V*theta_H,
            lambda : k_f[36]*theta_W,
            lambda : k_f[37]*theta_F*free,
            lambda : k_b[38]*0,                         #do-nothing formula for #38.
            lambda : k_f[39]*theta_I*free*free,
            lambda : k_f[40]*theta_B*free*free,
            lambda : k_f[41]*theta_CH3O*theta_H,
            lambda : k_f[42]*theta_CH*theta_H,
            lambda : k_f[43]*theta_CH2*theta_H,
            lambda : k_f[44]*theta_CH3*theta_H,
            lambda : k_f[45]*theta_OH*theta_H,
            lambda : k_f[46]*theta_CHO*theta_H,
            lambda : k_f[47]*theta_CH2O*theta_H,
            lambda : k_f[48]*theta_S,                   #phenol desorption
            lambda : k_f[49]*theta_M,                   #catechol desorption
            lambda : k_f[50]*theta_X,                   #benzene desorption
            lambda : k_f[51]*theta_CH3OH,               #CH3OH desorption
            lambda : k_f[52]*theta_CH4,                 #CH4 desorption
            lambda : k_f[53]*theta_H2O,                 #H2O desorption
            lambda : k_f[54]*theta_K,                   #anisole desorption
            #Note that in the original code, #55 below was listed as #71.
            #This was due to the numbering scheme favored by the chemist in his diagrams.
            #We changed the numbering scheme to better fit our programming needs here, 
            #thus indexes 55-61 correspond to indexes 71-77 in the original C++ code.
            lambda : k_f[55]*theta_O*free,
            lambda : k_f[56]*theta_Z*free,
            lambda : k_f[57]*theta_Z*free,
            lambda : k_f[58]*theta_Z1,
            lambda : k_f[59]*theta_G*theta_H,
            lambda : k_f[60]*theta_Y*free,
            lambda : k_f[61]*theta_CO*theta_H
        ]
    
        reverseRates = [
            lambda : k_b[0]*theta_A,
            lambda : k_b[1]*theta_B*free,
            lambda : k_b[2] *theta_C*free,          #special
            lambda : k_b[3] *theta_D*theta_OH,      #special
            lambda : k_b[4]*theta_E*theta_CH3O,
            lambda : k_b[5]*theta_F*theta_H,
            lambda : k_b[6]*theta_G*theta_CH3,
            lambda : k_b[7]*theta_I*theta_H,
            lambda : k_b[8]*theta_J*theta_H,
            lambda : k_b[9] *theta_K*theta_OH,      #special
            lambda : k_b[10]*theta_K*free,          #special
            lambda : k_b[11]*theta_S*free,
            lambda : k_b[12]*theta_E*theta_CH2O,
            lambda : k_b[13]*theta_L*theta_H,
            lambda : k_b[14]*theta_G*theta_CH2,
            lambda : k_b[15]*theta_M*free,
            lambda : k_b[16]*theta_N*theta_OH,
            lambda : k_b[17]*theta_N*theta_CH3O,
            lambda : k_b[18]*theta_O*theta_H,
            lambda : k_b[19]*theta_P*theta_CH2,
            lambda : k_b[20]*theta_Q*theta_H,       #special
            lambda : k_b[21]*theta_E*theta_CHO,
            lambda : k_b[22]*theta_G*theta_CH,
            lambda : k_b[23]*theta_E*theta_OH,
            lambda : k_b[24]*theta_T,
            lambda : k_b[25]*theta_R*free,
            lambda : k_b[26]*theta_N*theta_CH2O,
            lambda : k_b[27]*theta_U*theta_CH2,
            lambda : k_b[28]*theta_T*free,
            lambda : k_b[29]*theta_R*theta_CH2,     #special
            lambda : k_b[30]*theta_S*free,
            lambda : k_b[31]*theta_V*theta_OH,
            lambda : k_b[32]*theta_W*free,
            lambda : k_b[33]*theta_S*theta_OH,
            lambda : k_b[34]*theta_G*free*free,
            lambda : k_b[35]*theta_X*free*free,
            lambda : k_b[36]*theta_X*theta_OH,
            lambda : k_b[37]*theta_O*theta_H,
            lambda : k_b[38]*0,                         #do-nothing formula for #38.
            lambda : k_b[39]*theta_U*theta_CH3,
            lambda : k_b[40]*theta_S*theta_CH3O,
            lambda : k_b[41]*theta_CH3OH*free,
            lambda : k_b[42]*theta_CH2*free,
            lambda : k_b[43]*theta_CH3*free,
            lambda : k_b[44]*theta_CH4*free,
            lambda : k_b[45]*theta_H2O*free,
            lambda : k_b[46]*theta_CH2O*free,
            lambda : k_b[47]*theta_CH3O*free*free,
            lambda : k_b[48]*p_S*free*free*free*free,   #phenol desorption
            lambda : k_b[49]*p_M*free*free*free*free,   #catechol desorption
            lambda : k_b[50]*p_X*free*free*free,        #benzene desorption
            lambda : k_b[51]*p_CH3OH*free,              #CH3OH desorption
            lambda : k_b[52]*p_CH4*free,                #CH4 desorption
            lambda : k_b[53]*p_H2O*free,                #H2O desorption
            lambda : k_b[54]*p_K*free*free*free*free,   #anisole desorption
            #Note that in the original code, #55 below was listed as #71.
            #This was due to the numbering scheme favored by the chemist in his diagrams.
            #We changed the numbering scheme to better fit our programming needs here, 
            #thus indexes 55-61 correspond to indexes 71-77 in the original C++ code.
            lambda : k_b[55]*theta_Z*theta_H,
            lambda : k_b[56]*theta_U*theta_CH,
            lambda : k_b[57]*theta_Z1*theta_H,
            lambda : k_b[58]*theta_N*theta_CO,
            lambda : k_b[59]*theta_Y*free,
            lambda : k_b[60]*theta_R*theta_OH,
            lambda : k_b[61]*theta_CHO
        ]
    
        #Total rates of our reactions, forward_rate - reverse_rate.
        #This list of functions is generated by evaluating the forward and reverse rates
        #and then putting those values into a basic "a - b" sort of function.
        r = [(lambda a, b: a() - b())(f, r) for (f, r) in zip(forwardRates, reverseRates)]


        #Set global variables 'fr', 'rr', and 'tr' to the results from our forward, reverse, and total rates, respectively.
        global fr
        global rr
        global tr
        
        fr = [f() for f in forwardRates]
        rr = [f() for f in reverseRates]
        tr = r[:]

        
        #Note that in the original code, #55 below was listed as #71.
        #This was due to the numbering scheme favored by the chemist in his diagrams.
        #We changed the numbering scheme to better fit our programming needs here, 
        #thus indexes 55-61 correspond to indexes 71-77 in the original C++ code.

        #The steady state approximation should give us the resulting coverages.
        #These functions apply the steady state approximation for all thetas. This is our output to the solver.
        #When the output values of all of these functions are 0, the solver knows it's found the roots.
        steadyStateApprox = [
	        lambda r: r[4] + r[12] + r[21] + r[23] - r[11],     # d(theta_E)/dt=0
	        lambda r: r[5] - r[12] - r[13] - r[14] - r[37],     # d(theta_F)/dt=0
	        lambda r: r[6] + r[14] + r[22] + r[34] - r[15] - r[16] - r[59], # d(theta_G)/dt=0
	        lambda r: r[7] - r[17] - r[18] - r[39],            # d(theta_I)/dt=0
	        lambda r: r[13] - r[21] - r[22],                   # d(theta_L)/dt=0
	        lambda r: r[15] - r[23] - r[24] - r[49],            # d(theta_M)/dt=0 catechol
	        lambda r: r[16] + r[17] + r[26] + r[58] - r[25],     # d(theta_N)/dt=0
	        lambda r: r[18] + r[37] - r[26] - r[27] - r[55],     # d(theta_O)/dt=0
	        lambda r: r[25] + r[29] + r[60] - r[30],            # d(theta_R)/dt=0
	        lambda r: r[11] + r[30] + r[33] + r[40] - r[31] - r[32] - r[48], # d(theta_S)/dt=0 phenol
	        lambda r: r[27] + r[39] + r[56] - r[34],            # d(theta_U)/dt=0
	        lambda r: r[31] - r[35],                          # d(theta_V)/dt=0
	        lambda r: r[35] + r[36] - r[50],                   # d(theta_X)/dt=0 benzene

	        lambda r: r[3] + r[9] + r[23] + r[31] + r[31] + r[33] - r[45], # d(theta_OH)/dt=0
	        lambda r: r[45] - r[53],                          # d(theta_H2O)/dt=0
	        lambda r: r[22] + r[56] - r[42],                   # d(theta_CH)/dt =0
	        lambda r: r[14] + r[19] + r[27] + r[29] + r[42] - r[43], # d(theta_CH2)/dt=0	
	        lambda r: r[6] + r[43] - r[44],                   # d(theta_CH3)/dt=0
	        lambda r: r[44] - r[52],                          # d(theta_CH4)/dt=0
	        lambda r: r[21] + r[61] - r[46],                   # d(theta_CHO)/dt=0
	        lambda r: r[12] + r[26] + r[46] - r[47],            # d(theta_CH2O)/dt=0
	        lambda r: r[4] + r[17] + r[47] - r[41],            # d(theta_CH3O)/dt=0
	        lambda r: r[41] - r[51],                          # d(theta_CH3OH)/dt=0
	        lambda r: r[0] - r[1] - r[4] - r[5] - r[6] - r[7], # d(theta_A)/dt=0
        #	lambda r: theta_A - k_00_f/k_00_b*p_A*free*free*free*free; # d(theta_A)/dt=0 if step 0 is in equilibrium.
	        lambda r: r[32] - r[36],                          # d(theta_W)/dt=0
	        lambda r: r[24] + r[28] - r[33],                   # d(theta_T)/dt=0
	        lambda r: r[1] - r[8] - r[40],                   # d(theta_B)/dt=0
	        lambda r: r[8] - r[19],                          # d(theta_J)/dt=0
	        lambda r: r[19] - r[28],                          # d(theta_P)/dt=0
	        lambda r: theta_H - exp(-( (-1.374+0.076+0.683)/2 + 2*0.084*(theta_H-0.139))/KB/T)*pow(p_H2,0.5)*free, #theta_H
	        lambda r: r[2] - r[9],                          # d(theta_C)/dt=0
	        lambda r: r[3] - r[10],                          # d(theta_D)/dt=0
	        lambda r: r[9] + r[10] - r[20] - r[54],            # d(theta_K)/dt=0
	        lambda r: r[20] - r[29],                          # d(theta_Q)/dt=0
	        lambda r: r[59] - r[60],                          # d(theta_Y)/dt=0
	        lambda r: r[55] - r[56] - r[57],                   # d(theta_Z)/dt=0
	        lambda r: r[57] - r[58],                          # d(theta_Z1)/dt=0
	        lambda r: theta_CO - exp(-(-2.131+0.028+1.764)/KB/T)*p_CO*free, # d(theta_CO)/dt=0
                lambda r: 4*theta_A+4*theta_E+4*theta_F+4*theta_G+4*theta_I+4*theta_L+4*theta_M+
                4*theta_N+4*theta_O+4*theta_R+4*theta_S+5*theta_U+3*theta_V+3*theta_X+free+theta_H+
                theta_OH+theta_H2O+theta_CH+2*theta_CH2+theta_CH3+theta_CH4+2*theta_CHO+2*theta_CH2O+
                theta_CH3O+theta_CH3OH+4*theta_W+5*theta_T+4*theta_B+4*theta_J+5*theta_P+4*theta_C+
                4*theta_D+4*theta_K+5*theta_Q+4*theta_Y+4*theta_Z+4*theta_Z1+theta_CO-1.00
        ]


        #Set global variable 'ss' to the values calculated in our steady-state approximation:
        global ss
        ss = [f(r) for f in steadyStateApprox]

        #The end of our function:
        return [f(r) for f in steadyStateApprox]



    #These are the constraints on the system, they will include both bounds constraints and conservation constraints.
    constraints = [
        lambda theta_A, theta_E, theta_F, theta_G, theta_I, theta_L, theta_M, theta_N, theta_O, theta_R,
               theta_S, theta_U, theta_V, theta_X, free, theta_H, theta_OH, theta_H2O, theta_CH, theta_CH2,
               theta_CH3, theta_CH4, theta_CHO, theta_CH2O, theta_CH3O, theta_CH3OH, theta_W, theta_T, theta_B, theta_J,
               theta_P, theta_C, theta_D, theta_K, theta_Q, theta_Y, theta_Z, theta_Z1, theta_CO: (minBounds[0] <= theta_A <= maxBounds[0]) == True,
        lambda theta_A, theta_E, theta_F, theta_G, theta_I, theta_L, theta_M, theta_N, theta_O, theta_R,
               theta_S, theta_U, theta_V, theta_X, free, theta_H, theta_OH, theta_H2O, theta_CH, theta_CH2,
               theta_CH3, theta_CH4, theta_CHO, theta_CH2O, theta_CH3O, theta_CH3OH, theta_W, theta_T, theta_B, theta_J,
               theta_P, theta_C, theta_D, theta_K, theta_Q, theta_Y, theta_Z, theta_Z1, theta_CO: (minBounds[1] <= theta_E <= maxBounds[1]) == True,
        lambda theta_A, theta_E, theta_F, theta_G, theta_I, theta_L, theta_M, theta_N, theta_O, theta_R,
               theta_S, theta_U, theta_V, theta_X, free, theta_H, theta_OH, theta_H2O, theta_CH, theta_CH2,
               theta_CH3, theta_CH4, theta_CHO, theta_CH2O, theta_CH3O, theta_CH3OH, theta_W, theta_T, theta_B, theta_J,
               theta_P, theta_C, theta_D, theta_K, theta_Q, theta_Y, theta_Z, theta_Z1, theta_CO: (minBounds[2] <= theta_F <= maxBounds[2]) == True
    ]

    #Actually solve our system of equations:
    #NLSystem = findroot([coverage] + constraints, initialState, solver='secant', verbose=True, verify=False)
    NLSystem = findroot(coverage, initialState, solver='mnewton', verbose=False, verify=False)


    #--------------------------------------
    # DISPLAY RESULTS:

    i = 0
    print("FORWARD REACTION RATES:")
    for x in fr:
        print("r_f["+str(i)+"] = "+nstr(x))
        i += 1

    i = 0
    print("REVERSE REACTION RATES:")
    for x in rr:
        print("r_b["+str(i)+"] = "+nstr(x))
        i += 1

    i = 0
    print("TOTAL REACTION RATES:")
    for x in tr:
        print("r["+str(i)+"] = "+nstr(x))
        i += 1

    thetasList = [
        'theta_A', 'theta_E', 'theta_F', 'theta_G', 'theta_I', 'theta_L', 'theta_M', 'theta_N', 'theta_O', 'theta_R',
        'theta_S', 'theta_U', 'theta_V', 'theta_X', 'free', 'theta_H', 'theta_OH', 'theta_H2O', 'theta_CH', 'theta_CH2',
        'theta_CH3', 'theta_CH4', 'theta_CHO', 'theta_CH2O', 'theta_CH3O', 'theta_CH3OH', 'theta_W', 'theta_T', 'theta_B', 'theta_J',
        'theta_P', 'theta_C', 'theta_D', 'theta_K', 'theta_Q', 'theta_Y', 'theta_Z', 'theta_Z1', 'theta_CO'
    ]

    i = 0
    print("THETAS:")
    for x in NLSystem:
        print(thetasList[i]+" = "+nstr(x))
        i += 1

    i = 0
    print("STEADY STATE APPROX.:")
    for x in ss:
        print("out["+str(i)+"] = "+nstr(x))
        i += 1















