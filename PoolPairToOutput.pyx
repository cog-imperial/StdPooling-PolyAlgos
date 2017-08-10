import itertools

class PoolO(object):
    @staticmethod
    def PoolPairAtOneOutputSol(int o, pools,allOutputsConcUps,allOutputFlowUps,allOutputProfits,allInputsConc, allInputCosts,
                               inputsForPools,directsForOutputs, int TOL):
        # if direct coincides with input to second pool, choose direct - less coupling
        cdef double P_U = allOutputsConcUps[o];
        cdef double P_L = 0;
        cdef double D_U = round(allOutputFlowUps[o] * TOL, 0);
        cdef double d   = round(allOutputProfits[o] * TOL, 0);
        ######### find all 2 pool solutions and their profits analytically
        cdef int pool1 = pools[0];
        cdef int pool2 = pools[1];
        inputsToPool1 = list(set(inputsForPools[pool1]).difference(inputsForPools[pool2]+directsForOutputs[o]));
        inputsToPool2 = list(set(inputsForPools[pool2]).difference(inputsForPools[pool1]+directsForOutputs[o]));
        inputsToPool1 = list(zip(inputsToPool1,
                                  [allInputsConc[i] for i in inputsToPool1],
                                  [allInputCosts[i] for i in inputsToPool1]));
        inputsToPool2 = list(zip(inputsToPool2,
                                  [allInputsConc[i] for i in inputsToPool2],
                                  [allInputCosts[i] for i in inputsToPool2]));
        inputPairsList = list(itertools.product([elem for elem in inputsToPool1 if elem[1] < P_U],
                                          [elem for elem in inputsToPool2 if elem[1] > P_U])) + \
                   list(itertools.product([elem for elem in inputsToPool1 if elem[1] > P_U],
                                          [elem for elem in inputsToPool2 if elem[1] < P_U])) + \
                   list(itertools.product([elem for elem in inputsToPool1 if elem[1] < P_L],
                                          [elem for elem in inputsToPool2 if elem[1] > P_L])) + \
                   list(itertools.product([elem for elem in inputsToPool1 if elem[1] > P_L],
                                          [elem for elem in inputsToPool2 if elem[1] < P_L]));
        fSols = [0]*len(inputPairsList);
        for ind, (i, j) in enumerate(inputPairsList):
            P = P_L if ((i[1] - j[1]) * (i[2] - j[2]) > 0 and (i[1] - P_L) * (j[1] - P_L) < 0) or (
                (P_L <= j[1] <= P_U) and (i[1] <= P_L)) else P_U;
            xi = D_U * (round(P * TOL, 0) - round(j[1] * TOL, 0)) / ( round(i[1] * TOL, 0) - round(j[1] * TOL, 0));
            xj = D_U * (round(P * TOL, 0) - round(i[1] * TOL, 0)) / ( -round(i[1] * TOL, 0) + round(j[1] * TOL, 0));
            f = (d * D_U - round(i[2]*TOL,0) * xi - round(j[2]*TOL,0) * xj) / (TOL ** 2);
            #fSols[ind] = [pools, f, [i[1],j[1]], [i[0]], [xi/TOL], [j[0]], [xj/TOL] ];
            fSols[ind] = f;
        if not fSols:
            return 0;
        return max(fSols);