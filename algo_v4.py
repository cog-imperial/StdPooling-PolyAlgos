import math
from numpy import isreal
import numpy as np
from scipy import poly1d, roots
from sympy import symbols, Poly
from gams import GamsWorkspace,GamsException,GamsExceptionExecution
import os,glob
import sys
from collections import Counter
import cython_I_1_J_unconstr as cy
from PoolPairToOutput import PoolO
#import PoolPairToOutput
import operator
import itertools
from functools import reduce
import pyximport; pyximport.install()
from pyomo.environ import *
import pyomo.opt.solver
import time, copy

try:
    from line_profiler import LineProfiler


    def do_profile(follow=[]):
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()

            return profiled_func

        return inner

except ImportError:
    def do_profile(follow=[]):
        "Helpful if you accidentally leave in production!"

        def inner(func):
            def nothing(*args, **kwargs):
                return func(*args, **kwargs)

            return nothing

        return inner

class I_1_J_unconstr:
    def __init__(self, Cx, gx, TZ, Cz, gz, P_L, P_U, D_U, d, TX=[], tol=100, dec=6):
        self.TZ = TZ
        self.Cx = [round(x * tol, 0) for x in Cx];
        self.Cz = [round(x * tol, 0) for x in Cz]
        self.gx = [round(x * tol, 0) for x in gx]
        self.gz = [round(x * tol, 0) for x in gz]
        self.P_L = [round(x * tol, 0) for x in P_L]
        self.P_U = [round(x * tol, 0) for x in P_U]
        self.D_U = [round(x * tol, 0) for x in D_U]
        self.d = [round(x * tol, 0) for x in d]
        self.TX = TX
        self.tol = tol
        self.dec = dec

    def solve_I_1_1(self, o=0):
        tol = self.tol;
        dec = self.dec;
        Cx = self.Cx;
        gx = self.gx;
        Cz = self.Cz;
        gz = self.gz;  # Cz = [self.Cz[index] for index in TZ]; gz = [self.gz[index] for index in TZ];
        P_L = self.P_L[o];
        P_U = self.P_U[o];
        d = self.d[o];
        D_U = self.D_U[o];
        TX = self.TX;
        Z_sol = self.find_Z_sol(o);
        XZ_active = self.find_XZ_active(o);
        XZ_sol = [];  # store objective value, breakpoint, non-zero
        bps = [XZ_active[0][1]] + [x[2] for x in XZ_active];  # get a unique list of breakpoints
        for bp in bps:
            # get the first active set around the breakpoint
            active_set = [x[0] for x in XZ_active if x[1] == bp or x[2] == bp];
            active_set = active_set[0];
            if active_set == 'infeasible': continue;  # if input-only active set is infeasible and cannot be made feasible with any direct
            if len(active_set) == 2:  # if input-only active set at breakpoint
                i = active_set[0];
                j = active_set[1];
                xi = D_U * (bp - Cx[j]) / (Cx[i] - Cx[j]);
                xj = D_U * (bp - Cx[i]) / (Cx[j] - Cx[i]);
                f = (d * D_U - gx[i] * xi - gx[j] * xj) / (tol ** 2);
                if xj == 0:
                    XZ_sol.append([f, bp / tol, [TX[i]], [xi / tol], [], []]);
                elif xi == 0:
                    XZ_sol.append([f, bp / tol, [TX[j]], [xj / tol], [], []]);
                else:
                    XZ_sol.append([f, bp / tol, [TX[i], TX[j]], [xi / tol, xj / tol], [], []]);
            else:  # if mixed active set at breakpoint
                i = active_set[0];
                j = active_set[1];
                q = active_set[2];
                gXPair = (gx[i] * (bp - Cx[j]) + gx[j] * (Cx[i] - bp)) / (Cx[i] - Cx[j]);
                Pq = P_L if ((bp - Cz[q]) * (gXPair - gz[q]) > 0 and (bp - P_L) * (Cz[q] - P_L) < 0) or (
                    (P_L <= Cz[q] <= P_U) and (bp <= P_L)) else P_U;
                xi = D_U * (bp - Cx[j]) * (Pq - Cz[q]) / ((Cx[i] - Cx[j]) * (bp - Cz[q]));
                xj = D_U * (bp - Cx[i]) * (Pq - Cz[q]) / ((Cx[j] - Cx[i]) * (bp - Cz[q]));
                zq = D_U * (bp - Pq) / (bp - Cz[q]);
                f = (d * D_U - gx[i] * xi - gx[j] * xj - gz[q] * zq) / (tol ** 2);
                zIndices = [];
                zFlows = [];
                if xi < 0 or xj < 0:
                    print('Negative flow mixed pair!!!');
                if zq != 0: zIndices = [q]; zFlows = [zq / tol];
                if xj == 0:
                    XZ_sol.append([f, bp / tol, [TX[i]], [xi / tol], zIndices, zFlows]);
                elif xi == 0:
                    XZ_sol.append([f, bp / tol, [TX[j]], [xj / tol], zIndices, zFlows]);
                else:
                    XZ_sol.append(
                        [f, bp / tol, [TX[i], TX[j]], [xi / tol, xj / tol], zIndices, zFlows]);
        XZ_sol = XZ_sol + [Z_sol];
        maxProfit = max(x[0] for x in XZ_sol);
        return [x for x in XZ_sol if x[0] == maxProfit];  # XZ_sol;#

    def solve_I_1_J(self):
        tol = self.tol;
        dec = self.dec;
        Cz = self.Cz;
        activeForAllO = [];
        bpsCommon = [];
        for o in range(0, len(self.TZ)):
            # [Z_sol, Z_Cost] = self.find_Z_sol(o,1);
            XZ_active = self.find_XZ_active(o,
                                            1);  # get all active sets for output o together with their p-derivatives
            bps = [XZ_active[0][1]] + [x[2] for x in XZ_active];  # get a unique list of breakpoints for output o
            activeForAllO.append(XZ_active);
            bpsCommon = bpsCommon + bps;
        # bpsCommon = list(unique_everseen(bpsCommon)); # get unique list of bps
        bpsCommon = sorted(list(set(bpsCommon)));
        XZsolsAtBps = [];
        for i in range(0, len(bpsCommon) - 1):  # for all pairs of breakpoints common among all outputs
            bl = bpsCommon[i];
            bu = bpsCommon[i + 1];  # lower and upper breakpoints of the interval
            derivsMixed = [];  # list of derivatives for all outputs in the interval where a mixed set dominates
            derivTotalInput = 0;
            for o in range(0, len(self.TZ)):  # for each output
                [active_set, bp, bp2, deriv] = [XZ for XZ in activeForAllO[o] if XZ[1] <= bl and bu <= XZ[2]][0];
                # if direct solution or no feasible solution found, no derivative influencing p
                if active_set[0] == 'direct' or active_set[0] == 'infeasible': continue;
                if len(active_set) == 2:
                    derivTotalInput = derivTotalInput + deriv;
                else:
                    derivsMixed.append([deriv, Cz[
                        active_set[
                            2]] / tol]);  # append derivative factor plus Cq - still need to divide by (p-Cq)^2
            nbDerivsPos = sum(x[0] > 0 for x in derivsMixed) + (derivTotalInput > 0) * 1;
            # find solutions at interval endpoints
            XZsolsAtBps.append(
                self.findTotalFlowsAtConcInInterval(activeForAllO, bl, bl, bu));  # solution at bl concentration
            XZsolsAtBps.append(
                self.findTotalFlowsAtConcInInterval(activeForAllO, bu, bl, bu));  # solutiion at bu concentration
            # if all mixed derivs and the sum of all input-only derivs are all positive or all negative
            # it means the max is at one interval endpoint and there is no need to solve
            # if not, we need to solve univariate polynomial for p where stationary point/possibly max may occur
            if not (nbDerivsPos == 0 or nbDerivsPos == (len(derivsMixed) + (derivTotalInput != 0) * 1)):
                p = symbols('p', real=True);
                derivsMixed = [[deriv, (p - Cq) ** 2] for [deriv, Cq] in derivsMixed];
                mixedTerms = sum([deriv * self.productFactors(
                    pFactor for (idx2, [deriv2, pFactor]) in enumerate(derivsMixed) if idx2 != idx1) \
                                  for (idx1, [deriv, Cq]) in enumerate(derivsMixed)]);
                inputTerm = self.productFactors(pFactor for [deriv2, pFactor] in derivsMixed) * derivTotalInput;
                expr = poly1d(Poly(mixedTerms + inputTerm, p).all_coeffs());
                p_sols = [root * tol for root in roots(expr) if isreal(root) and bl < root * tol < bu];
                for p in p_sols: XZsolsAtBps.append(self.findTotalFlowsAtConcInInterval(activeForAllO, p, bl, bu));
        maxProfit = max(x[0] for x in XZsolsAtBps);
        return [x for x in XZsolsAtBps if x[0] == maxProfit] + [XZsolsAtBps];

    def find_Z_sol(self, o=0, costAlso=False):
        tol = self.tol;
        dec = self.dec;
        Z_sol = [];
        TZ = self.TZ[o];
        Cz = self.Cz;
        gz = self.gz;
        P_L = self.P_L[o];
        P_U = self.P_U[o];
        d = self.d[o];
        D_U = self.D_U[o];
        # find lowest cost direct node
        lowestCost, index = min((cost, index) for (index, cost) in enumerate(gz) if (index in TZ));
        if P_L <= Cz[index] <= P_U:  # if feasible
            if d - lowestCost > 0:  # if profitable
                Z_sol = [D_U * (d - lowestCost) / (tol ** 2), [], [], [], [index], [D_U / tol, dec]];
        else:
            # find set of all distinct feasible pairs of direct inputs dominating both their nodes
            s = [];
            for i in range(0, len(TZ)):
                for j in range(i + 1, len(TZ)):
                    gi = gz[TZ[i]];
                    gj = gz[TZ[j]];
                    Ci = Cz[TZ[i]];
                    Cj = Cz[TZ[j]];
                    if (gi < gj and not (P_L <= Ci <= P_U) and P_L <= Cj <= P_U) \
                            or (gi > gj and P_L <= Ci <= P_U and not (
                                            P_L <= Cj <= P_U)) or Ci < P_L < Cj or Ci > P_L > Cj:
                        Pij = P_L if (Ci - Cj) / (gi - gj) > 0 else P_U;
                        g = (gi * (Pij - Cj) + gj * (Ci - Pij)) / (Ci - Cj);
                        s.append([g, TZ[i], TZ[j], Pij]);
            if s:
                # find the dominant direct-only active set and its solution
                lowestCost, index = min((pair[0], index) for (index, pair) in enumerate(s));
                i = s[index][1];
                j = s[index][2];
                Pij = s[index][3];
                zi = D_U * (Pij - Cz[j]) / (Cz[i] - Cz[j]);
                zj = D_U * (Pij - Cz[i]) / (Cz[j] - Cz[i]);
                fz = d * D_U - zi * gz[i] - zj * gz[j];
                if zi < 0 or zj < 0:
                    print('Negative flow direct pair!!!');
                if zi == 0:
                    zIndices = [j];
                    zFlows = [zj / tol];
                elif zj == 0:
                    zIndices = [i];
                    zFlows = [zi / tol];
                else:
                    zIndices = [i, j];
                    zFlows = [zi / tol, zj / tol];
                Z_sol = [fz / (tol ** 2), [], [], [], zIndices, zFlows];
        if not Z_sol: return [[], []];  # no feasible direct-only solution available
        if costAlso:
            # directFlows = [0] * len(TZ); # for storing and then adding all direct flows across outputs
            # qsInTZ = [index for (index, TZo) in enumerate(TZ) if TZo in Z_sol[4]];
            # for i in range(0,len(qsInTZ)):
            #    directFlows[qsInTZ[i]] = Z_sol[5][i];
            return [[Z_sol[0], Z_sol[4], Z_sol[5]], lowestCost];
        else:
            return Z_sol;

    def find_XZ_active(self, o=0, derivs=False):
        tol = self.tol;
        dec = self.dec;
        Cx = self.Cx;
        gx = self.gx;
        TZ = self.TZ[o];
        Cz = self.Cz;
        gz = self.gz;
        P_L = self.P_L[o];
        P_U = self.P_U[o];
        D_U = self.D_U[o];
        d = self.d[o];
        X_bps = self.find_X_bps();
        X_bps = self.splitXintervalsAroundQualityBounds(X_bps, P_L, P_U);
        XZ_active = [];
        P2 = [];  # store in P2 all pairs of distinct direct pairs
        for i in range(0, len(TZ)):
            for j in range(i + 1, len(TZ)): P2.append([TZ[i], TZ[j]]);
        for X in X_bps:  # for each input pair in its breakpoint interval [Il,Iu]
            i = X[0][0];
            j = X[0][1];
            Il = X[1];
            Iu = X[2];
            gi = gx[i];
            gj = gx[j];
            Ci = Cx[i];
            Cj = Cx[j];
            # 1,2. find R and Q
            Q = [];
            R = [];
            gl = (gi * (Il - Cj) + gj * (Ci - Il)) / (Ci - Cj);
            gu = (gi * (Iu - Cj) + gj * (Ci - Iu)) / (Ci - Cj);
            if P_L <= Iu <= P_U and P_L <= Il <= P_U:
                Q = [z for z in TZ if (Cz[z] < P_L or Cz[z] > P_U) and gz[z] < max(gl, gu)];  # mixed triple 1
                if not Q:
                    XZ_active.append(X);
                    continue;  # input pair dominates entire interval
                else:
                    R = Q;
            else:
                if Iu > P_U:
                    Q = [z for z in TZ if Cz[z] < P_U];  # mixed triple 2
                else:
                    Q = [z for z in TZ if Cz[z] > P_L];  # mixed triple 3
                if not Q:
                    XZ_active.append(['infeasible', Il, Iu]);
                    continue;  # input/mixed sets infeasible
                else:
                    R = [q for q in Q if gz[q] < max(gl, gu)];
            P2feasible = [pair for pair in P2 if
                          (pair[0] in Q and pair[1] in Q)];  # find all pairs of feasible mixed sets;
            # 3. truncate interval X into S and B
            S = [Il, Iu];
            Bs = [];
            if gi == gj:  # in case of cost equality of the two inputs
                cheaperDirects = [q for q in R if gz[q] < gi];
                if cheaperDirects: S = [];
            else:
                slope = (gi - gj) / (Ci - Cj);
                for q in R:
                    b = (Ci * (gz[q] - gj) - Cj * (gz[q] - gi)) / (gi - gj);
                    if slope > 0:
                        if b > S[0]:
                            S[1] = b;
                        else:
                            S = [];
                            break;
                    else:
                        if b < S[1]:
                            S[0] = b;
                        else:
                            S = [];
                            break;
            if not S:
                Bs.append([Il, Iu]);  # at the end of truncating
            elif S == [Il, Iu]:
                Bs = [];
            elif S[0] == Il:
                Bs.append([S[1], Iu]);
            elif S[1] == Iu:
                Bs.append([Il, S[0]]);
            else:
                Bs.append([Il, S[0]]);
                Bs.append([S[1], Iu]);
            a1 = gi - gj;
            a2 = gi * Cj - gj * Ci;
            a3 = Ci - Cj;
            # 4. find dominant mixed sets on intervals B where input pair X is dominated by mixed triples
            if Bs:
                for B in Bs:
                    PB = self.getDominantMixedActiveSetsOnInterval(B, R, P2feasible, Il, [a1, a2, a3], 0, o);
                    # new_bps = self.getActiveSetsFromBreakpoints(B,PB,X[0]);
                    # new_bps = [bp P_L if xDominates*(Cq>Il) else P_U for bp in bps];
                    XZ_active = XZ_active + self.getActiveSetsFromBreakpoints(B, PB, X[0]);
            # 5. find dominant mixed sets on intervals S where input pair X dominates any mixed triples (but is feasible/infeasible)
            if S and P_L <= Iu <= P_U and P_L <= Il <= P_U:  # X feasible
                XZ_active.append([[i, j], S[0], S[1]]);
            elif S:  # X infeasible
                PS = self.getDominantMixedActiveSetsOnInterval(S, Q, P2feasible, Il, [a1, a2, a3], 1, o);
                XZ_active = XZ_active + self.getActiveSetsFromBreakpoints(S, PS, X[0]);

        # if we are solving a multiple outputs problem and need derivatives and breakpoints w.r.t. direct-only set also
        if derivs:
            XZ_active2 = [];
            for idx, [active_set, bp, bp2] in enumerate(XZ_active):  # calculate df(p)/dp for all active sets found
                #idx=idx+7;
                [Z_sol, gZ] = self.find_Z_sol(o, True);
                if active_set == 'infeasible':  # if input-only active set is infeasible but there is a feasible direct-only solution
                    if not Z_sol:
                        XZ_active2.append([['infeasible'], bp, bp2, []]);
                    else:
                        XZ_active2.append([['direct'] + Z_sol[1], bp, bp2, Z_sol]);
                    continue;
                deriv = 0;
                p = 0;
                i = active_set[0];
                j = active_set[1];
                Ci = Cx[i];
                Cj = Cx[j];
                gi = gx[i];
                gj = gx[j];
                # cost for input-only pair at bp (lower end of interval)
                gXZAtbp = (gi * (bp - Cj) + gj * (Ci - bp)) / (Ci - Cj);
                if len(active_set) == 2:  # if input-only active set between breakpoints
                    deriv = -D_U * (gi - gj) / (Ci - Cj) / tol;
                    if Z_sol:
                        if gi != gj:
                            # potential breakpoint between input-only pair and direct-only set
                            p = round((Ci * (gZ - gj) - Cj * (gZ - gi)) / (gi - gj), dec);
                        else:
                            p = 0;
                else:  # if mixed active set between breakpoints
                    q = active_set[2];
                    Cq = Cz[q];
                    gq = gz[q];
                    Pq = P_L if ((bp - Cq) * (gXZAtbp - gq) > 0 and (bp - P_L) * (Cq - P_L) < 0) or (
                        (P_L <= Cq <= P_U) and (bp2 <= P_L)) else P_U;
                    deriv = -D_U * (Pq - Cq) / (Ci - Cj) * (Cq * (gj - gi) + Cj * (gi - gq) + Ci * (gq - gj)) / (
                        tol ** 3);
                    if Z_sol:
                        # cost for mixed triple at bp (lower end of interval)
                        gXZAtbp = (gq * (Pq - bp) + gXZAtbp * (Cq - Pq)) / (Cq - bp);
                        # potential breakpoint between mixed triple and direct-only set
                        pNum = (gq * Pq - gZ * Cq) * (Ci - Cj) + (gj * Ci - gi * Cj) * (Cq - Pq);
                        pDenom = (gq - gZ) * (Ci - Cj) - (gi - gj) * (Cq - Pq);
                        # special case: if direct q in triple is the same as direct in Z_sol, if Z_sol is better choose it
                        if q == Z_sol[1][0] and pDenom==0:
                            if gXZAtbp < gZ:
                                XZ_active2.append([active_set, bp, bp2, deriv]); continue;
                            else:
                                XZ_active2.append([['direct'] + Z_sol[1], bp, bp2, Z_sol]); continue;
                        p = round(pNum / pDenom, dec);
                # if no feasible direct-only set exists, add existing active set and its deriv
                if not Z_sol: XZ_active2.append([active_set, bp, bp2, deriv]); continue;
                # check whether the potential breakpoint materializes in both input vs direct and mixed vs direct cases
                if (bp < p < bp2):  # if breakpoint p in interval
                    if gXZAtbp < gZ:
                        XZ_active2.append([active_set, bp, p, deriv]);  # adjust interval for input only pair
                        XZ_active2.append([['direct'] + Z_sol[1], p, bp2,
                                           Z_sol]);  # add also direct set on correct interval with deriv 0;
                    else:
                        XZ_active2.append([['direct'] + Z_sol[1], bp, p,
                                           Z_sol]);  # add also direct set on correct interval with deriv 0;
                        XZ_active2.append([active_set, p, bp2, deriv]);  # adjust interval for input only pair
                else:
                    if gXZAtbp < gZ:
                        XZ_active2.append([active_set, bp, bp2, deriv]);
                    else:
                        XZ_active2.append([['direct'] + Z_sol[1], bp, bp2, Z_sol]);
            # sort XZ_active_2 by first breakpoint
            XZ_active2.sort(key=lambda x: x[1]);
            # merge breakpoint intervals where directs dominate
            XZ_active3 = [];
            prevIdx = -10;
            directInterv = [];
            for idx, [active_set, bp, bp2, deriv] in enumerate(XZ_active2):
                if active_set[0] == 'direct':
                    if idx == prevIdx + 1:
                        directInterv = [active_set, directInterv[1], bp2, deriv];
                    else:
                        directInterv = [active_set, bp, bp2, deriv];
                    if idx < (len(XZ_active2) - 1):  # if not at the last elements
                        prevIdx = idx;  # see if one identical direct comes next
                    else:
                        XZ_active3.append(directInterv);  # if last element, append the direct-only set
                elif prevIdx != -10:
                    XZ_active3.append(directInterv);
                    XZ_active3.append([active_set, bp, bp2, deriv]);
                    prevIdx = -10;
                else:
                    XZ_active3.append([active_set, bp, bp2, deriv]);

            # find infeasible breakpoints and their intervals, where sending the flow to the output is not profitable
            # this will only happen for input-only or mixed sets, since we eliminated unprofitable sets in Z_sol directly (constant profit function)
            XZ_active4 = [];
            for idx, [active_set, bp, bp2, deriv] in enumerate(XZ_active3):
                if active_set[0] == 'direct' or active_set[0] == 'infeasible': XZ_active4.append(
                    [active_set, bp, bp2, deriv]); continue;
                i = active_set[0];
                j = active_set[1];
                Ci = Cx[i];
                Cj = Cx[j];
                gi = gx[i];
                gj = gx[j];
                isProfitableAtBp = -1;
                pNum = -1;
                pDenom = -1;
                if len(active_set) == 2:  # for input-only sets
                    # calculate potential profitability/feasibility breakpoint p (the numerator part)
                    pNum = d * (Ci - Cj) + gi * Cj - gj * Ci;
                    pDenom = gi - gj;
                    # calculate sign of profit at lower bp
                    isProfitableAtBp = (d - (gi * (bp - Cj) - gj * (bp - Ci)) / (Ci - Cj)) > 0;
                else:  # for mixed input pairs
                    q = active_set[2];
                    gXPair = (gx[i] * (p - Cx[j]) + gx[j] * (Cx[i] - p)) / (Cx[i] - Cx[j]);
                    Pq = P_L if ((p - Cz[q]) * (gXPair - gz[q]) > 0 and (p - P_L) * (Cz[q] - P_L) < 0) or (
                        (P_L <= Cz[q] <= P_U) and (bp2 <= P_L)) else P_U;
                    pNum = -d * (Ci - Cj) * Cq + (Pq - Cq) * (gi * Cj - gj * Ci);
                    pDenom = (gi - gj) * (Pq - Cq) - d * (Ci - Cj);
                    isProfitableAtBp = (
                                           d - (Pq - Cq) * (gi * (bp - Cj) - gj * (bp - Ci)) / (
                                           (Ci - Cj) * (bp - Cq))) > 0;
                # check profitability conditions for either input or mixed pair
                if bp * pDenom < pNum < bp2 * pDenom and pDenom != 0:  # if there is a valid profitability breakpoint
                    p = pNum / pDenom;
                    if isProfitableAtBp:
                        XZ_active4.append([active_set, bp, p, deriv]);
                        XZ_active4.append([['infeasible'], p, bp2, []]);
                    else:
                        XZ_active4.append([['infeasible'], bp, p, []]);
                        XZ_active4.append([active_set, p, bp2, deriv]);
                else:
                    if isProfitableAtBp:
                        XZ_active4.append([active_set, bp, bp2, deriv]);
                    else:
                        XZ_active4.append([['infeasible'], bp, bp2, []]);
                        # merge breakpoint intervals with infeasibility that are next to each other
            XZ_active5 = [];
            prevIdx = -10;
            infeasInterv = [];
            for idx, [active_set, bp, bp2, deriv] in enumerate(XZ_active4):
                if (active_set[0] == 'infeasible'):
                    if idx == prevIdx + 1:
                        infeasInterv = [active_set, infeasInterv[1], bp2, deriv];
                    else:
                        infeasInterv = [active_set, bp, bp2, deriv];
                    if idx < (len(XZ_active2) - 1):  # if not at the last elements
                        prevIdx = idx;  # see if one identical infeas comes next
                    else:
                        XZ_active5.append(infeasInterv);  # if last element, append the infeas set
                elif prevIdx != -10:
                    XZ_active5.append(infeasInterv);
                    XZ_active5.append([active_set, bp, bp2, deriv]);
                    prevIdx = -10;
                else:
                    XZ_active5.append([active_set, bp, bp2, deriv]);
            return XZ_active5;
        return XZ_active;

    def find_X_bps(self):
        X_bps = [];
        Cx = self.Cx;
        gx = self.gx;
        lC = len(Cx);
        P1 = [];  # store in P1 all pairs of distinct input pairs, their obj slope and intercept
        for i in range(0, lC):
            for j in range(i + 1, lC):
                P1.append((i, j));
        # P2 = np.reshape(np.transpose(np.meshgrid(np.arange(lC), np.arange(lC))), (lC * lC, 2));
        bps = [];  # find all breakpoints occuring at one input concentration
        for i in range(0, len(Cx)):
            isBp = True;
            p = Cx[i];
            # gxarr = np.array(gx);
            # Cxarr = np.array(Cx);
            # Cxarr0 = Cxarr[P2[:, 0]]; Cxarr1=Cxarr[P2[:, 1]];
            # gAgregNomArr = gxarr[P2[:, 0]] * (p - Cxarr1) + gxarr[P2[:, 1]] * (Cxarr0 - p);
            # gAgregDenomArr = Cxarr0 - Cxarr1;
            # isBp = np.logical_not(np.any((P2[:, 0] != i) * (P2[:, 1] != i) * ((p - Cxarr0) * (p - Cxarr1) < 0) * \
            #     (gAgregNomArr < gx[i] * gAgregDenomArr) * (gAgregDenomArr > 0)));
            for P in P1:
                p0 = P[0];
                p1 = P[1];
                if p0 != i and p1 != i and (p - Cx[p0]) * (
                    p - Cx[p1]) < 0:  # ((Cx[P[0]] <= p <= Cx[P[1]]) or (Cx[P[0]] >= p >= Cx[P[1]])):
                    # gAgreg = (gx[P[0]]*(p-Cx[P[1]]) + gx[P[1]]*(Cx[P[0]]-p)) / (Cx[P[0]]-Cx[P[1]]);
                    gAgregNom = gx[p0] * (p - Cx[p1]) + gx[p1] * (Cx[p0] - p);
                    gAgregDenom = Cx[p0] - Cx[p1];
                    if (gAgregNom < gx[i] * gAgregDenom) * (gAgregDenom > 0) or \
                                    (gAgregNom > gx[i] * gAgregDenom) * (gAgregDenom < 0): isBp = False; continue;
            if isBp: bps.append([i, p]);
        bps.sort(key=lambda x: x[1]);
        # put in format [active pair, bp1, bp2]
        for bp in range(0, len(bps) - 1):
            X_bps.append([[bps[bp][0], bps[bp + 1][0]], bps[bp][1], bps[bp + 1][1]]);
        return X_bps;

    def getDominantMixedActiveSetsOnInterval(self, B, R, P2, Il, a, xDominates, o=0):
        Cz = self.Cz;
        gz = self.gz;
        P_L = self.P_L[o];
        P_U = self.P_U[o];
        PB = [];
        bl = B[0];
        bu = B[1];
        a1 = a[0];
        a2 = a[1];
        a3 = a[2];
        for P in P2:  # for each pair of mixed sets {x,q}, {x,r}
            gq = gz[P[0]];
            gr = gz[P[1]];
            Cq = Cz[P[0]];
            Cr = Cz[P[1]];
            if gq == gr: continue;  # for equal costs of direct, one mixed set will dominate throughout
            Pq = P_L if (xDominates * (Cq > Il) and (Il - P_L) * (Cq - P_L) < 0) or (
                (P_L <= Cq <= P_U) and (bu <= P_L)) else P_U;
            Pr = P_L if (xDominates * (Cr > Il) and (Il - P_L) * (Cr - P_L) < 0) or (
                (P_L <= Cr <= P_U) and (bu <= P_L)) else P_U;
            # find 1 or 2 breakpoints between mixed sets {x,q}, {x,r}
            b1 = Pq - Cq - Pr + Cr;
            b2 = Pr * Cq - Pq * Cr;
            c1 = gq - gr;
            c2 = gr * (Cq + Pr) - gq * (Cr + Pq);
            c3 = gq * Pq * Cr - gr * Pr * Cq;
            ae = a1 * b1 + a3 * c1;
            be = a1 * b2 - a2 * b1 + a3 * c2;
            ce = -b2 * a2 + a3 * c3;
            delta = be ** 2 - 4 * ae * ce;
            if delta < 0: continue;
            p1 = (-be + math.sqrt(delta)) / (2 * ae);
            p2 = (-be - math.sqrt(delta)) / (2 * ae);
            for p in ([p1, p2] if delta else [p1]):  # for each potential breakpoint/root within B
                if bl < p < bu:
                    # gxp  = (gi*(p-Cj) +  gj*(Ci-p) ) / (Ci-Cj);
                    gxp = (-a2 + p * a1) / a3;
                    gxqp = (gq * (Pq - p) + gxp * (Cq - Pq)) / (Cq - p);
                    # if q and r direct nodes are cheaper than X pair, or X pair does not dominate
                    if gxp > max(gq, gr) or not xDominates:
                        # keep mixed pair if dominant above all other mixed pairs
                        mixedCosts = [];
                        for w in R:
                            Pw = P_L if (xDominates * (Cz[w] > Il) and (Il - P_L) * (Cz[w] - P_L) < 0) or (
                                (P_L <= Cz[w] <= P_U) and (bu <= P_L)) else P_U;
                            mixedCosts.append((gz[w] * (Pw - p) + gxp * (Cz[w] - Pw)) / (Cz[w] - p));
                        if gxqp == min(mixedCosts):
                            PB.append([p, P[0], P[1]]);
        gxbl = (-a2 + bl * a1) / a3;
        mixedCosts = [];
        for w in R:
            Pw = P_L if (xDominates * (Cz[w] > Il) and (Il - P_L) * (Cz[w] - P_L) < 0) or (
                (P_L <= Cz[w] <= P_U) and (bu <= P_L)) else P_U;
            mixedCosts.append([(gz[w] * (Pw - bl) + gxbl * (Cz[w] - Pw)) / (Cz[w] - bl), w]);
        gxqbl, q = min((cost, index) for (cost, index) in mixedCosts);
        if not PB:
            PB.append([q]);  # if no breakpoint occurs add dominating mixed pair
        else:
            PB.sort(key=lambda x: x[0]);  # sort PB by breakpoint concentration p
            if q == PB[0][2]:  # for the lowest/first breakpoint order the mixed sets
                PB[0] = [PB[0][0], PB[0][2], PB[0][1]];
        return PB;

    def findTotalFlowsAtConcInInterval(self, activeForAllO, p, bl, bu):
        tol = self.tol;
        dec = self.dec;
        Cx = self.Cx;
        gx = self.gx;
        Cz = self.Cz;
        gz = self.gz;
        totalInputFlows = [0] * len(Cx);
        totalDirectFlows = [];
        totalObj = 0;
        TX = self.TX;
        totalDirectFlows = [[[], []]] * len(self.TZ);  # for storing and then adding all direct flows across outputs
        for o in range(0, len(self.TZ)):  # for each output
            [active_set, _, _, Z_sol] = [XZ for XZ in activeForAllO[o] if XZ[1] <= bl and bu <= XZ[2]][0];
            if active_set[0] == 'infeasible':
                continue;
            elif active_set[0] == 'direct':
                # return [[Z_sol[0], directFlows, Z_sol[4]], lowestCost];
                # return [Z_sol[0], Z_sol[4], Z_sol[5], lowestCost];
                totalDirectFlows[o] = Z_sol[1:3];
                f = Z_sol[0];
            else:
                P_L = self.P_L[o];
                P_U = self.P_U[o];
                d = self.d[o];
                D_U = self.D_U[o];
                xi = 0;
                xj = 0;
                i = active_set[0];
                j = active_set[1];
                if len(active_set) == 2:
                    xi = D_U * (p - Cx[j]) / (Cx[i] - Cx[j]);
                    xj = D_U * (p - Cx[i]) / (Cx[j] - Cx[i]);
                    f = round((d * D_U - gx[i] * xi - gx[j] * xj) / (tol ** 2), dec);
                else:
                    q = active_set[2];
                    gXPair = (gx[i] * (p - Cx[j]) + gx[j] * (Cx[i] - p)) / (Cx[i] - Cx[j]);
                    Pq = P_L if ((p - Cz[q]) * (gXPair - gz[q]) > 0 and (p - P_L) * (Cz[q] - P_L) < 0) or (
                        (P_L <= Cz[q] <= P_U) and (bu <= P_L)) else P_U;
                    xi = D_U * (p - Cx[j]) * (Pq - Cz[q]) / ((Cx[i] - Cx[j]) * (p - Cz[q]));
                    xj = D_U * (p - Cx[i]) * (Pq - Cz[q]) / ((Cx[j] - Cx[i]) * (p - Cz[q]));
                    zq = D_U * (p - Pq) / (p - Cz[q]);
                    f = round((d * D_U - gx[i] * xi - gx[j] * xj - gz[q] * zq) / (tol ** 2), dec);
                    totalDirectFlows[o] = [[q], [round(zq / tol, dec)]];
                totalInputFlows[i] = totalInputFlows[i] if xi == 0 else totalInputFlows[i] + round(xi / tol, dec);
                totalInputFlows[j] = totalInputFlows[j] if xj == 0 else totalInputFlows[j] + round(xj / tol, dec);
            totalObj = totalObj + f;
        return [totalObj, p / tol, [TX[index] for index, flow in enumerate(totalInputFlows) if flow != 0],
                [flow for flow in totalInputFlows if flow != 0], totalDirectFlows];

    @staticmethod
    def splitXintervalsAroundQualityBounds(X_bps, P_L, P_U):
        for X in X_bps:
            removeX = False;
            if X[1] < P_L < X[2]:
                removeX = True;
                X_bps.append([X[0], X[1], P_L]);
                X_bps.append([X[0], P_L, X[2]]);
            if X[1] < P_U < X[2]:
                removeX = True;
                X_bps.append([X[0], X[1], P_U]);
                X_bps.append([X[0], P_U, X[2]]);
            if removeX: X_bps.remove(X);
        return X_bps;

    @staticmethod
    def getActiveSetsFromBreakpoints(I, PI, Xpair):
        pl = I[0];
        pu = I[1];
        # If no bp on I, a pair dominates with bps at I bounds, if 1 bp get AI directly
        if not PI: return [[Xpair, pl, pu]];
        P0 = PI[0];
        if len(PI) == 1:
            if len(P0) == 1: return [
                [Xpair + P0, pl, pu]];  # if no breakpoint occurs add only dominating mixed pair
            return [[Xpair + [P0[1]], pl, P0[0]], [Xpair + [P0[2]], P0[0], pu]];
        AI = [[Xpair + [P0[1]], pl, P0[0]]];  # add active set on first (lowest) bp interval
        # PI ordered w.r.t. bps; At each bp (start from 2nd) order active sets by previous bp
        for i in range(1, len(PI)):
            Pi = PI[i];
            Pi_prev = PI[i - 1];
            if Pi_prev[2] != Pi[1]:
                AI.append([Xpair + [Pi_prev[2]], Pi_prev[0], Pi[0]])  # add active set on intermediate bp interval
                Pi[2] = Pi[1];  # partial swap, as only Pi[2] is used at next iteration
            if i == len(PI) - 1:
                AI.append([Xpair + [Pi[2]], Pi[0], pu]);  # add active set on last (highest) bp interval
        return AI;

    @staticmethod
    def productFactors(list):
        r = 1;
        for x in list: r *= x;
        return r;

# if __name__ == "__main__":
#timeinOld=0;
def formatAndRun2PoolsGams(o, lines3, ws, probCostMatrix, pools, inputNb,poolNb,outputNb,
                           allPoolFlowUps,allOutputFlowUps, allOutputsConcUps, allInputsConcNp, solPairsForO):
    for idx, line in enumerate(lines3):
        if "c(i,j)" in line:
            lines3[idx + 1] = str('    \t    %d\t   %d\t    %d\n' % (inputNb + 1, inputNb + 2, inputNb + 3));
            lines3[idx + 2] = "\n".join(("\t".join((str("%.2f" % val).rjust(6) if ij > 0 else str("%d" % val).rjust(3)
                                                    for ij, val in enumerate(x))) for x in probCostMatrix)) + ';\n';
            del lines3[idx + 3: idx + 2 + inputNb + poolNb];
        if "a(i,j)" in line:
            lines3[idx + 1] = str('    \t %d\t %d\t %d\n' % (inputNb + 1, inputNb + 2, inputNb + 3));
            # probAdjMatrix = probCostMatrix.copy();
            probCostMatrix[:, 1:][probCostMatrix[:, 1:] != 0] = 1;
            lines3[idx + 2] = "\n".join(
                ("\t".join((str('%d' % _).rjust(3) for _ in x)) for x in probCostMatrix)) + ';\n';
            del lines3[idx + 3: idx + 2 + inputNb + poolNb];
        if "q(i,k)" in line:  # or
            lines3[idx + 2 + inputNb] = str(inputNb + 3).rjust(3) + '\t' + str(allOutputsConcUps[o - 1]).rjust(
                5) + ';\n';
            del lines3[idx + 3 + inputNb: idx + 2 + inputNb + outputNb];
        if "bl(i)" in line:
            lines3[idx + inputNb] = str('\t\t  %d\t%.2f' % (inputNb + 3, 0)) + ' /;\n';
            del lines3[idx + 1 + inputNb: idx + inputNb + outputNb];
        if "bu(i)" in line:
            lines3[idx + inputNb:idx + inputNb + 1] = \
                [str('\t\t  %d\t%.2f\n' % (inputNb + 1 + ij, allPoolFlowUps[pool])) for ij, pool in enumerate(pools)];
            lines3[idx + inputNb + len(pools)] = \
                str('\t\t  %d\t%.2f /;\n' % (inputNb + 3, allOutputFlowUps[o - 1]));
            del lines3[idx + inputNb + 1 + len(pools): -1];
    allLines = "".join(lines3);
    t2 = ws.add_job_from_string(allLines);
    success = 0;
    while not success:
        try:
            t2.run();
            success = 1;
        except GamsExceptionExecution:
            pass
    poolsAndDirectsActive = [int(rec.key(0)) - 1 for rec in t2.out_db["f"] if rec.level != 0];
    if (inputNb in poolsAndDirectsActive) and (inputNb + 1 in poolsAndDirectsActive):
        # if both pools are active
        # poolsAndDirectsActive.remove(inputNb); poolsAndDirectsActive.remove(inputNb+1);
        inputsActive = [int(rec.key(0)) - 1 for rec in t2.out_db["y"] if rec.level != 0];
        inputsActive.sort();
        lower = 0;
        upper = 0;
        for rec in t2.out_db.get_symbol("profit"):
            lower = int(rec.lower);
            upper = int(rec.upper);
        inputsDirectsForEachPool = [sorted(s[3] + s[5]) for s in solPairsForO if s[0] in pools];
        maxProfitEachPool = max([s[1] for s in solPairsForO if s[0] in pools] + [0]);
        if inputsActive not in inputsDirectsForEachPool and maxProfitEachPool < upper:
            p1 = sum(allInputsConcNp[[int(rec.key(0)) - 1 for rec in t2.out_db["y"] if
                                      rec.level != 0 and int(rec.key(1)) - 1 == inputNb]] *
                     [rec.level for rec in t2.out_db["y"] if rec.level != 0 and int(rec.key(1)) - 1 == inputNb]);
            p2 = sum(allInputsConcNp[[int(rec.key(0)) - 1 for rec in t2.out_db["y"] if
                                      rec.level != 0 and int(rec.key(1)) - 1 == inputNb + 1]] *
                     [rec.level for rec in t2.out_db["y"] if rec.level != 0 and int(rec.key(1)) - 1 == inputNb + 1]);
            inputs1 = [int(rec.key(0)) - 1 for rec in t2.out_db["x"] if
                       int(rec.key(1)) - 1 == inputNb and int(rec.level) != 0];
            flows1 = [rec.level for rec in t2.out_db["x"] if int(rec.key(1)) - 1 == inputNb and int(rec.level) != 0];
            inputs2 = [int(rec.key(0)) - 1 for rec in t2.out_db["x"] if
                       int(rec.key(1)) - 1 == inputNb + 1 and int(rec.level) != 0];
            flows2 = [rec.level for rec in t2.out_db["x"] if
                      int(rec.key(1)) - 1 == inputNb + 1 and int(rec.level) != 0];
            return [pools, lower, [p1, p2], inputs1, flows1, inputs2, flows2];
        return [];

def formatAndRunGams(lines3, ws, probCostMatrix, pools, outputs, inputNb,poolNb,outputNb,
                           allOutputFlowUps, allOutputsConcUps):
    pools = list(pools);
    for idx, line in enumerate(lines3):
        # replace sets cardinalities
        if "set i" in line:
            lines3[idx] = '    set i    / 1*{:d}/;\n'.format(inputNb + len(pools)+len(outputs));
        if "set t(i)" in line:
            lines3[idx] = '    set t(i) / {:d}*{:d}/;\n'.format(inputNb + len(pools)+1, inputNb + len(pools)+len(outputs));
        if "c(i,j)" in line:
            lines3[idx + 1] = str(('\t   %d'*(len(pools)+len(outputs)) % tuple(range(inputNb+1,inputNb+len(pools)+len(outputs)+1))) + "\n");
            lines3[idx + 2] = "\n".join(("\t".join((str("%.2f" % val).rjust(6) if ij > 0 else str("%d" % val).rjust(3)
                                                    for ij, val in enumerate(x))) for x in probCostMatrix)) + ';\n';
            del lines3[idx + 3: idx + 2 + inputNb + poolNb];
        if "a(i,j)" in line:
            lines3[idx + 1] = str(('\t %d' * (len(pools)+len(outputs))
                                         % tuple(range(inputNb + 1, inputNb + len(pools)+ len(outputs) + 1))) + '\n');
            probCostMatrix[:, 1:][probCostMatrix[:, 1:] != 0] = 1;
            lines3[idx + 2] = "\n".join(
                ("\t".join((str('%d' % _).rjust(3) for _ in x)) for x in probCostMatrix)) + ';\n';
            del lines3[idx + 3: idx + 2 + inputNb + poolNb];
        if "q(i,k)" in line:  # or
            for ix,output in enumerate(outputs):
                lines3[idx + 2 + inputNb+ix] = str(inputNb + len(pools)+ix+1).rjust(3) + '\t' + \
                                               str(allOutputsConcUps[output]).rjust(5) + '\n';
                if ix == len(outputs)-1:
                    lines3[idx + 2 + inputNb + ix] = str(inputNb + len(pools) + ix + 1).rjust(3) + '\t' + \
                                                     str(allOutputsConcUps[output]).rjust(5) + ';\n';
            del lines3[idx + 2 + inputNb+len(outputs): idx + 2 + inputNb + outputNb];
        if "bl(i)" in line:
            for ix in range(0,len(pools)+len(outputs)+1):
                lines3[idx + inputNb+ix] = str('\t\t  %d\t%.2f' % (inputNb + ix+1, 0)) + '\n';
                if ix == len(pools)+len(outputs) - 1:
                    lines3[idx + inputNb + ix] = str('\t\t  %d\t%.2f' % (inputNb + ix + 1, 0)) + ' /;\n';
            del lines3[idx + len(pools+outputs)+ inputNb: idx + inputNb + outputNb];
        if "bu(i)" in line:
            lines3[idx + inputNb:idx + inputNb + len(pools)] = \
                    [str('\t\t  %d\t+inf\n' % (inputNb + 1 + ij)) for ij, pool in enumerate(pools)];
            lines3[idx + inputNb + len(pools):idx + inputNb + len(pools)+len(outputs)-1] = \
                [str('\t\t  %d\t%.2f \n' % (inputNb + 1+len(pools)+ij, allOutputFlowUps[output])) for ij, output in enumerate(outputs) if ij<len(outputs)-1];
            lines3[idx + inputNb + len(pools) + len(outputs)-1] = \
                str('\t\t  %d\t%.2f /;\n' % (inputNb + len(pools) + len(outputs), allOutputFlowUps[outputs[len(outputs)-1]]) );
            del lines3[idx + inputNb +len(pools)+len(outputs): -1];
    allLines = "".join(lines3);
    t2 = ws.add_job_from_string(allLines);
    success = 0;
    #tt0 = time.time();
    while not success:
        try:
            t2.run();
            for rec in t2.out_db.get_symbol("profit"):
                lower = float(rec.lower);
            if lower<10000000:
                success = 1;
        except GamsExceptionExecution:
            pass
    #tt1 = time.time();
    #global timeinOld
    #timeinOld += tt1 - tt0;
    genTime=0; solveTime=0; nodesUsed=0;
    for rec in t2.out_db.get_symbol("profit"):
        lower = float(rec.lower);
    for rec in t2.out_db.get_symbol("genTime"):
        genTime = rec.value;
    for rec in t2.out_db.get_symbol("solveTime"):
        solveTime = rec.value;
    for rec in t2.out_db.get_symbol("nodesUsed"):
        nodesUsed = rec.value;
    return [lower,genTime,solveTime,nodesUsed];

def formatAndRunBARON(lines2, ws):
    allLines = "".join(lines2);
    t0 = time.time();
    t2 = ws.add_job_from_string(allLines);
    success = 0; lowerB = 0; upperB =0;
    while not success:
        try:
            t2.run();
            for rec in t2.out_db.get_symbol("profit"):
                lowerB = float(rec.lower);
                upperB = float(rec.upper);
            if lowerB < 10000000:
                success = 1;
        except GamsExceptionExecution:
            pass
    t1 = time.time();
    timeBARON = t1-t0;
    genTime = 0;
    solveTime = 0;
    nodesUsed = 0;
    for rec in t2.out_db.get_symbol("profit"):
        lower = float(rec.lower);
    for rec in t2.out_db.get_symbol("genTime"):
        genTime = rec.value;
    for rec in t2.out_db.get_symbol("solveTime"):
        solveTime = rec.value;
    for rec in t2.out_db.get_symbol("nodesUsed"):
        nodesUsed = rec.value;
    return [lowerB,upperB,timeBARON,genTime,solveTime,nodesUsed];

def PoolPairAtOneOutputSol2(o,pools,allOutputsConcUps,allOutputFlowUps,allOutputProfits,allInputsConc, allInputCosts,
                           inputsForPools,directsForOutputs, TOL):
            # if direct coincides with input to second pool, choose direct - less coupling
            P_U = allOutputsConcUps[o];
            P_L = 0;
            D_U = round(allOutputFlowUps[o] * TOL, 0);
            d   = round(allOutputProfits[o] * TOL, 0);
            ######### find all 2 pool solutions and their profits analytically
            pool1 = pools[0];
            pool2 = pools[1];
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
            return max(fSols);#max(fSols, key=lambda x: x[1]);

def findPoolGroups(poolpairs):
    #poolpairs = [[2, 20], [2, 3], [4, 2], [3, 5], [4, 5], [4, 12], [12, 13], [7, 14], [14, 22], [100, 101]];
    coupledLists = [];
    while poolpairs:
        coupledPairs = poolpairs[0];
        coupledList = poolpairs[0];
        while coupledPairs:
            coupledPairs = [p[1] if p[0] in coupledList else p[0] for p in poolpairs
                            if (p[0] in coupledList and p[1] not in coupledList) or (
                            p[1] in coupledList and p[0] not in coupledList)];
            coupledList = coupledList + list(set(coupledPairs));
            poolpairs = [p for p in poolpairs if not (p[0] in coupledList and p[1] in coupledList)];
        coupledLists.append(tuple(coupledList));
    return coupledLists;

def assocPoolsConcs(poolsOutputsAssign):
    poolsAndSolConcs=[];
    for elem in poolsOutputsAssign:
        if len(elem[0]) == 1:
            if elem[0][0] != -1:
                poolsAndSolConcs.append((elem[0][0], elem[2]));
        else:
            poolsAndSolConcs.append((elem[0][0], elem[2][0]));
            poolsAndSolConcs.append((elem[0][1], elem[2][1]));
    return list(set(poolsAndSolConcs));


#@do_profile(follow=[])
def main():
    t0=time.time();
    TOL = 100;  # tolerance for floating point arithmetic precision
    DEC = 6;  # round to how many decimals to cut floating point errors
    problemName = "sppC3.gms";
    #problemName = "Haverly1.gms";
    # problemName = "Foulds3.gms";

    filelist = glob.glob("I:\OneDrive for Business 1\PhD Radu\Ruth-MILP\Python\Pooling Python Code\data\_gams_py*");
    for file in filelist:
        os.remove(file);

    inputNb = 0;
    outputNb = 0;
    poolNb = 0;
    kNb = 0;
    inputsForPools = [];
    directsForOutputs = [];
    poolsForOutputs = [];
    allInputCosts = [];
    allOutputProfits = [];
    allInputsKs = [];
    allOutputsKsUps = [];
    allPoolFlowUps = [];
    allOutputFlowUps = [];

    # parse problem from gams file
    # *************************************************************************
    lines = [];
    with open(os.path.join(os.path.curdir, "data/", problemName), "r") as dataFile:
        lines = dataFile.readlines();
    lines[4]="$offlisting";
    for i, line in enumerate(lines):
        if "set s(i)" in line:
            inputNb = int(line[line.index('*') + 1:line.index('/;')].strip());
        elif "set t(i)" in line:
            poolNb = int(line[line.index('/') + 1:line.index('*')].strip()) - 1 - inputNb;
            outputNb = int(line[line.index('*') + 1:line.index('/;')].strip()) - inputNb - poolNb;
        elif "set k" in line:
            kNb = int(line[line.index('*') + 1:line.index('/;') - 1].strip());
            lines[i] = "    set k /1*1/ ;\n";
            break;
    inputCostRows = [];
    poolCostRows = [];
    adjRows = [];
    for idx, line in enumerate(lines):
        if "c(i,j)" in line:
            inputCostRows = lines[idx + 2: idx + 2 + inputNb];
            poolCostRows = lines[idx + 2 + inputNb: idx + 2 + inputNb + poolNb];
        if "a(i,j)" in line:
            adjRows = lines[idx + 2: idx + 2 + inputNb + poolNb];
            break;
    inputCostMatrix = np.zeros((inputNb, poolNb + outputNb));
    inputOrigCostMatrix = np.zeros((inputNb, poolNb + outputNb));
    poolCostMatrix = np.zeros((poolNb, poolNb + outputNb));
    adjMatrix = np.zeros((inputNb + poolNb, poolNb + outputNb + 1));
    for idx, line in enumerate(poolCostRows):  poolCostMatrix[idx] = np.fromstring(line, dtype=float, sep=' ')[1:];
    for idx, line in enumerate(adjRows):       adjMatrix[idx] = np.fromstring(line, dtype=float, sep=' ');
    del adjRows;
    allOutputProfits = np.ndarray.tolist(-np.min(poolCostMatrix[:, poolNb:], axis=0));
    for idx, line in enumerate(inputCostRows):
        inputCostsPlusProfits = np.fromstring(line, dtype=float, sep=' ')[1:];
        bitMask = np.copy(inputCostsPlusProfits);
        bitMask[bitMask != 0] = 1;
        inputCostMatrix[idx] = inputCostsPlusProfits + bitMask * np.hstack((np.zeros(poolNb), allOutputProfits));
        inputOrigCostMatrix[idx] = inputCostsPlusProfits;
    allInputCosts = np.ndarray.tolist(np.max(inputCostMatrix, axis=1));
    for pool in range(0, poolNb): inputsForPools.append(np.nonzero(inputCostMatrix[:, pool])[0].tolist());
    for output in range(poolNb, poolNb + outputNb):
        directs = np.nonzero(inputCostMatrix[:, output])[0];
        directsForOutputs.append(directs.tolist());
        pools = np.nonzero(poolCostMatrix[:, output])[0];
        poolsForOutputs.append(pools.tolist());
    for idx, line in enumerate(lines):
        if "q(i,k)" in line:
            inputKRows = lines[idx + 2: idx + 2 + inputNb];
            outputKRows = lines[idx + 2 + inputNb: idx + 2 + inputNb + outputNb];
            allInputsKs = np.zeros((inputNb, kNb));
            allOutputsKsUps = np.zeros((outputNb, kNb));
            for idx, line in enumerate(inputKRows):      allInputsKs[idx] = np.fromstring(line, dtype=float, sep=' ')[1:];
            for idx, line in enumerate(outputKRows): allOutputsKsUps[idx] = np.fromstring(line, dtype=float, sep=' ')[1:];
        if "bu(i)" in line:
            poolFlowUpRows = lines[idx + inputNb : idx + inputNb + poolNb];
            for line in poolFlowUpRows:
                allPoolFlowUps.append(np.fromstring(line, dtype=float, sep=' ')[1]);
            outputFlowUpRows = lines[idx + inputNb + poolNb: idx + inputNb + poolNb + outputNb];
            for line in outputFlowUpRows:
                allOutputFlowUps.append(np.fromstring(line, dtype=float, sep=' ')[1]);
    for idx, line in enumerate(poolCostRows):  poolCostMatrix[idx] = np.fromstring(line, dtype=float, sep=' ')[1:];
    # ************************************************************end of parsing
    for idx, lineText in enumerate(lines):
        if "bu(i)" in lineText:
            InputPoolFlowUpRows = lines[idx: idx + inputNb + poolNb];
            for i, line in enumerate(InputPoolFlowUpRows):
                lines[idx + i] = line[0:line.rfind(' ')] + " +inf \n";

    if len(sys.argv) > 1:
        ws = GamsWorkspace(system_directory=sys.argv[1])
    else:
        ws = GamsWorkspace(working_directory="..\\Pooling Python Code\\data",debug=0) #ws = GamsWorkspace()
        #ws = GamsWorkspace(debug=DebugLevel.KeepFiles)
    t1 = time.time();

    filename = os.path.join(os.path.curdir, "data\\", "%s_results3.txt" % problemName);
    string = 'Reading in/interpreting file - time took = %f s' % (t1-t0);
    with open(filename, "a") as text_file: text_file.write(string+"\n"); print(string);
    # solve a separate problem for each quality k
    for k in range(2, kNb):  # kNb):
        ##### clean temporary gams files
        filelist = glob.glob(
            "I:\OneDrive for Business 1\PhD Radu\Ruth-MILP\Python\Pooling Python Code\data\_gams_py*");
        for file in filelist:
            os.remove(file);
        string = '*SOLVING problem with quality # %d' % (k+1);
        with open(filename, "a") as text_file: text_file.write('\n\n\n' + string+'\n'); print('\n\n'+string);
        allInputsConcNp = allInputsKs[:, k];
        allInputsConc = allInputsConcNp.tolist();
        allOutputsConcUps = allOutputsKsUps[:, k].tolist();
        lines2 = lines[:];

        # ************solve original problem in gams (modify to have one quality and no inputs/pools constraints)
        for idx, lineText in enumerate(lines2):
            if "q(i,k)" in lineText:
                lines2[idx + 1] = "\t    1\n";
                allConcs = allInputsConc + allOutputsConcUps;
                allIndices = list(range(1, inputNb + 1)) + list(
                    range(inputNb + poolNb + 1, inputNb + poolNb + outputNb + 1));
                for i, l in enumerate(allConcs):
                    lines2[idx + 2 + i] = str(allIndices[i]).rjust(3) + "\t" + str('%.2f' % l).rjust(5) + "\n";
                    lines2[idx + 1 + inputNb + outputNb] = str(allIndices[-1]) + "   " + str(
                    allConcs[inputNb + outputNb - 1]) + ";\n";
                del allConcs,allIndices;

        ###### SOLVE with BARON
        res = [0,0,0]; #formatAndRunBARON(lines2, ws)
        lowerB = res[0]; upperB =res[1]; timeBARON=res[2];

        string = 'BARON solve results: lower = %.6f ,  upper = %.6f - time= %f s' % (lowerB, upperB, timeBARON);
        with open(filename, "a") as text_file:
            text_file.write(string + '\n'); print(string);

        startTime=time.time();
        pairsForAllO = [];
        solPairsForAllO = [];
        upperBound = 0;
        poolsOutputsAssign = []; poolsAndSolConcs=[];
        ################## Create hierarchy for each output
        for o in range(1, outputNb+1):  # range(1,outputNb+1)
            #print("Creating hierarchy for output %d" % o);
            oColumn = adjMatrix[inputNb:inputNb + poolNb, poolNb + o];
            poolForO = [];
            poolPairsForO = [];
            for i in range(0, len(oColumn)):
                if oColumn[i]:
                    # add pair with only one pool/direct
                    poolForO.append(i);
                    # add pair of two pools/directs
                    for j in range(i + 1, len(oColumn)):
                        if oColumn[j]:
                            poolPairsForO.append([i, j]);
            pairsForAllO.append(poolForO);
            solPairsForO = [];
            # list of single pools that were inferior to directs (also 2 pool combination of them will be inferior)
            poolsWorseThanDirects = [];
            for pool in poolForO:
                #inputPairs = findOptimalInputPairs(allInputCosts,allInputsConc);
                #if len(pools) == 1:
                inputsToPool = inputsForPools[pool];
                prob = cy.I_1_J_unconstr( \
                    # inputs array to pool 1 with its concentrations and costs
                    [allInputsConc[i] for i in inputsToPool], [allInputCosts[i] for i in inputsToPool], \
                    # directs arrays to each output, followed by arrays of all directs concentrations and costs
                    [directsForOutputs[o - 1]], allInputsConc, allInputCosts, \
                    # lower and upper bound concetrations at each output
                    [0], [allOutputsConcUps[o - 1]], \
                    # upper limit of flow at each output
                    [allOutputFlowUps[o - 1]], \
                    # profit per unit flow sent to each output
                    [allOutputProfits[o - 1]], \
                    # inputs array to pool 1 with its indices
                    inputsToPool, TOL, DEC);
                ans = prob.solve_I_1_1();
                # prefer a direct solution if equal with one involving a pool, otherwise #TODO take all solutions
                if len(ans) > 1:
                    directsSol = [sol for sol in ans if sol[1] == []];
                    if directsSol:
                        solPairsForO.append([[-1]]+directsSol[0]);
                        poolsWorseThanDirects.append(pool);
                    else:
                        solPairsForO.append([[pool]]+ans[0]);
                else:
                    if ans[0][1] == []:
                        solPairsForO.append([[-1]]+ans[0]);
                        poolsWorseThanDirects.append(pool);
                    else:
                        solPairsForO.append([[pool]]+ans[0]);
            # traverse again in order to compare with 1 pool active sets
            # if direct coincides with input to second pool, choose direct - less coupling
            P_U = allOutputsConcUps[o - 1];
            P_L = 0;
            D_U = round(allOutputFlowUps[o - 1] * TOL, 0);
            d   = round(allOutputProfits[o - 1] * TOL, 0);
            solPoolPairForO=[];
            ############ UNCOMMENT code for finding sol for 2 pools, 1 output
            for pools in poolPairsForO:
                #pools=[1,8];
                if not set(pools).issubset(poolsWorseThanDirects) :
                    ######### find all 2 pool solutions and their profits analytically
                    pool1 = pools[0];
                    pool2 = pools[1];
                    inputsToPool1 = list(set(inputsForPools[pool1]).difference(inputsForPools[pool2]+directsForOutputs[o - 1]));
                    inputsToPool2 = list(set(inputsForPools[pool2]).difference(inputsForPools[pool1]+directsForOutputs[o - 1]));
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
                        fSols[ind] = [pools, f, [i[1],j[1]], [i[0]], [xi/TOL], [j[0]], [xj/TOL] ];
                    if not fSols:
                        continue;
                    fMax = max(fSols, key=lambda x: x[1]);
                    solMax = max(solPairsForO[poolForO.index(pool2)][1], solPairsForO[poolForO.index(pool1)][1]);
                    if fMax[1] - solMax> 0.000001:
                        solPoolPairForO.append(fMax);
            # remove duplicate entries with directs (-1)
            seen = set();
            solPairsForO = [item for item in solPairsForO if item[0][0] not in seen and not seen.add(item[0][0])];# remove duplicate entries with directs (-1)
            solPairsForO = solPairsForO+ solPoolPairForO; #omit for now treating 2 pools with multiple outputs (need to go to GAMS for that)
            solPairsForO.sort(key=lambda x: x[1], reverse=True);
            # if no pool provides flow to output at the best solution, choose this direct and separable solution
            if solPairsForO[0][0] == -1: solPairsForO = [solPairsForO[0]];
            solPairsForAllO.append(solPairsForO);
            upperBound = upperBound + solPairsForO[0][1];
            poolsOutputsAssign.append(solPairsForO[0]);
            if len(solPairsForO[0][0])==1:
                if solPairsForO[0][0][0]!=-1:
                    poolsAndSolConcs.append((solPairsForO[0][0][0],solPairsForO[0][2]));
            else:
                poolsAndSolConcs.append((solPairsForO[0][0][0], solPairsForO[0][2][0]));
                poolsAndSolConcs.append((solPairsForO[0][0][1], solPairsForO[0][2][1]));
            a=1;
        poolsAndSolConcs = list(set(poolsAndSolConcs));# unique
        t1 = time.time();
        string = 'Built output hierarchies - time= %f s' % (t1 - startTime);
        with open(filename, "a") as text_file:
            text_file.write(string+'\n'); print(string);
        t0 = time.time();
        iterCount = 0;
        seen=[];
        uniquePairPools = [item[0] for item in poolsOutputsAssign if item[0] != [-1] and
                           len(item[0]) == 2 and item[0] not in seen and not seen.append(item[0])];
        coupledPoolLists = findPoolGroups(uniquePairPools);
        seen = [item for sublist in coupledPoolLists for item in sublist];
        uniqueOnePools = [item[0][0] for item in poolsOutputsAssign if item[0]!=[-1] and
                          len(item[0]) == 1 and item[0][0] not in seen and not seen.append(item[0][0])];
        outputGroupSols = [0]*len(uniqueOnePools);
        # solve polynomially all one pool independent groups
        for i, pool in enumerate(uniqueOnePools):
            #### Solve multiple outputs problem in all original groups
            inputsToPool = inputsForPools[pool];
            outputs = [];
            for (output, item) in enumerate(poolsOutputsAssign):
                if item[0] == [pool]: outputs.append(output);
            # solve one pool, multiple outputs problem
            probOrig = cy.I_1_J_unconstr( \
                # inputs array to pool 1 with its concentrations and costs
                [allInputsConc[i] for i in inputsToPool], [allInputCosts[i] for i in inputsToPool], \
                # directs arrays to each output, followed by arrays of all directs concentrations and costs
                [directsForOutputs[o] for o in outputs], allInputsConc, allInputCosts, \
                # lower and upper bound concentrations at each output
                [0] * len(outputs), [allOutputsConcUps[o] for o in outputs], \
                # upper limit of flow at each output
                [allOutputFlowUps[o] for o in outputs], \
                # profit per unit flow sent to each output
                [allOutputProfits[o] for o in outputs], \
                # inputs array to pool 1 with its indices
                inputsToPool, TOL, DEC);
            outputGroupSols[i] = (pool, probOrig.solve_I_1_J()[0]);
        # solve using gams for any other multi-pool group and their associated outputs
        for pools in coupledPoolLists:
            if len(set(x for x in poolsAndSolConcs if x[0] in pools))==len(pools):
                outputGroupSols.append((pools,[sum([x[1] for x in poolsOutputsAssign if x[0][0] in pools])]));
            else:
                outputs= [ix for ix,item in enumerate(poolsOutputsAssign) if (item[0][0] in pools or (len(item[0])==2 and item[0][1] in pools))]
                reducedPoolCostMatrix = np.zeros(shape=(len(pools), len(outputs)));
                for ip,pool in enumerate(pools):
                    for io, output in enumerate(outputs):
                        poolAtOut = poolsOutputsAssign[output][0];
                        if pool==poolAtOut  or (type(poolAtOut) is list and pool in poolAtOut):
                            reducedPoolCostMatrix[ip,io] = poolCostMatrix[pool, poolNb + output];
                probCostMatrix = np.hstack(
                    [np.array([np.arange(1, len(inputOrigCostMatrix[:, 0]) + len(pools) + 1)]).T,
                     np.vstack([inputOrigCostMatrix[:, pools], np.zeros(shape=(len(pools), len(pools)))]),
                     np.vstack([inputOrigCostMatrix[:,poolNb+np.array(outputs)],reducedPoolCostMatrix]) ]);
                lines3 = lines2[:];
                profitPools = formatAndRunGams(lines3, ws, probCostMatrix, pools, outputs,
                    inputNb,poolNb,outputNb,allOutputFlowUps, allOutputsConcUps);
                outputGroupSols.append((pools,[profitPools[0]]));
            # CODE for only 2 POOLS AND 1 OUTPUT
            # for pools in poolPairsForO:
            #     # pools=[1,8];
            #     if not set(pools).issubset(poolsWorseThanDirects):
            #         probCostMatrix = np.hstack(
            #             [np.array([np.arange(1, len(inputOrigCostMatrix[:, 0]) + len(pools) + 1)]).T,
            #              np.vstack([inputOrigCostMatrix[:, pools], np.zeros(shape=(len(pools), len(pools)))]),
            #              np.array([np.hstack([np.zeros(inputNb),  # inputOrigCostMatrix[:, poolNb + o - 1],
            #                                   poolCostMatrix[pools, poolNb + o - 1]])]).T]);
            #         lines3 = lines2[:];
            #         sol1 = formatAndRun2PoolsGams(o, lines3, ws, probCostMatrix, pools,
            #                                       inputNb, poolNb, outputNb, allPoolFlowUps, allOutputFlowUps,
            #                                       allOutputsConcUps,
            #                                       allInputsConcNp, solPairsForO);
            #         sol2 = [];

        ######################### Iteration 0 !!!!!!
        objVal = sum([outputGroupSol[1][0] for outputGroupSol in outputGroupSols])\
                 +sum([x[1] for x in poolsOutputsAssign if x[0]==[-1]]);
        t1 = time.time();
        string ="Iteration %d, obj = %.6f\t - time= %f" % (iterCount, objVal, t1-t0)
        with open(filename, "a") as text_file:
            text_file.write(string+'\n');print(string);

        ##### loop in next iterations
        positiveMove = True; onlyCostsForAllO_prev=[]; poolsChanged=[]; savedSolves=0;
        while positiveMove: # while there are positive moves go through hierarchies for all outputs and make succesive moves
            movesForAllO = []; bestMovePerO=[0]*outputNb; onlyCostsForAllO=[];
            t0 = time.time();
            for (outp,solPairForO) in enumerate(solPairsForAllO): #for each output compute moves available and their cost/profit
                outp=outp+6;solPairForO=solPairsForAllO[outp];
                # if a direct dominates for output and no moves will occur here (profit of move is 0)
                if len(solPairForO)==1 and solPairForO[0][0] == -1: # when a direct dominates for output outp and there's no move
                    movesForAllO.append([[0,-1,-1]]); bestMovePerO[outp] = [0,-1,-1];
                    onlyCostsForAllO.append(0);
                    continue;
                if outp==17:
                    a=1;
                ######### Initial stage before move of pools choice for output
                initTotal = sum([s[1][0] for s in outputGroupSols]);

                ######### After dettaching current pools for output from the rest of the structure
                intTotal =0;
                poolOld = poolsOutputsAssign[outp][0]; # current pools being dettached;
                if type(poolOld) is not list: poolOld=[poolOld];
                poolsOutputsAssignInt = copy.deepcopy(poolsOutputsAssign);
                poolsOutputsAssignInt[outp][0]=[-1];
                seen = [];
                uniquePairPoolsInt = [item[0] for item in poolsOutputsAssignInt if item[0] != [-1] and
                                   len(item[0])==2 and item[0] not in seen and not seen.append(item[0])];
                coupledPoolListsInt = findPoolGroups(uniquePairPoolsInt);
                seen = [item for sublist in coupledPoolListsInt for item in sublist];
                uniqueOnePoolsInt = [item[0][0] for item in poolsOutputsAssignInt if item[0] != [-1] and
                                     len(item[0]) == 1 and item[0][0] not in seen and not seen.append(item[0][0])];
                outputGroupSolsInt=[];
                # add changed coupled pool groups
                for pools in [cp for cp in coupledPoolListsInt if set(cp)& set(poolOld)]: # for coupled pools sharing any pool with pooldOld pair
                    outputs = [ix for ix, item in enumerate(poolsOutputsAssignInt) if
                               (item[0][0] in pools or (len(item[0])==2 and item[0][1] in pools))];
                    revenue=0;
                    if len(outputs) == 1 and len(pools)==2:
                        revenue = PoolO.PoolPairAtOneOutputSol(
                            outputs[0], pools, allOutputsConcUps, allOutputFlowUps, allOutputProfits,
                            allInputsConc, allInputCosts, inputsForPools, directsForOutputs, TOL);
                    else:
                        reducedPoolCostMatrix = np.zeros(shape=(len(pools), len(outputs)));
                        for ip, pool in enumerate(pools):
                            for io, output in enumerate(outputs):
                                poolAtOut = poolsOutputsAssignInt[output][0];
                                if pool == poolAtOut or (type(poolAtOut) is list and pool in poolAtOut):
                                    reducedPoolCostMatrix[ip, io] = poolCostMatrix[pool, poolNb + output];
                        probCostMatrix = np.hstack(
                            [np.array([np.arange(1, len(inputOrigCostMatrix[:, 0]) + len(pools) + 1)]).T,
                             np.vstack([inputOrigCostMatrix[:, pools], np.zeros(shape=(len(pools), len(pools)))]),
                             np.vstack([inputOrigCostMatrix[:, poolNb + np.array(outputs)], reducedPoolCostMatrix])]);
                        lines3 = lines2[:];
                        res = formatAndRunGams(lines3, ws, probCostMatrix, pools, outputs,
                                        inputNb, poolNb, outputNb, allOutputFlowUps, allOutputsConcUps);
                        revenue = res[0];
                    intTotal += revenue;
                    outputGroupSolsInt.append((pools,[revenue]));
                # add changed one pool groups
                for pool in (set(uniqueOnePoolsInt) - set(uniqueOnePools)) | (set(uniqueOnePoolsInt) & set(poolOld)): # for sets of one pool that have been added or that share the pool with the poolOld
                    outputs=[];
                    for (output, item) in enumerate(poolsOutputsAssignInt):
                        if item[0] == [pool]: outputs.append(output);
                    #### Solve multiple outputs problem in the original group with the output outp removed and independent
                    inputsToPool = inputsForPools[pool];
                    probOrig = cy.I_1_J_unconstr( \
                        # inputs array to pool 1 with its concentrations and costs
                        [allInputsConc[i] for i in inputsToPool], [allInputCosts[i] for i in inputsToPool], \
                        # directs arrays to each output, followed by arrays of all directs concentrations and costs
                        [directsForOutputs[o] for o in outputs], allInputsConc, allInputCosts, \
                        # lower and upper bound concentrations at each output
                        [0] * len(outputs), [allOutputsConcUps[o] for o in outputs], \
                        # upper limit of flow at each output
                        [allOutputFlowUps[o] for o in outputs], \
                        # profit per unit flow sent to each output
                        [allOutputProfits[o] for o in outputs], \
                        # inputs array to pool 1 with its indices
                        inputsToPool, TOL, DEC);
                    revenue= probOrig.solve_I_1_J()[0][0];
                    intTotal += revenue;
                    outputGroupSolsInt.append((pool, [revenue]));
                changed = [cp for cp in coupledPoolListsInt if set(cp)& set(poolOld)];
                changed = [item for sublist in changed for item in sublist] +\
                        list((set(uniqueOnePoolsInt) - set(uniqueOnePools)) | (set(uniqueOnePoolsInt) & set(poolOld)));
                # add unchanged coupled or single pool groups
                for sol in outputGroupSols:
                    if type(sol[0]) is not tuple: pools=[sol[0]];
                    else:pools=sol[0];
                    if not set(pools) & set(changed) and not set(pools) & set(poolOld):
                        intTotal += sol[1][0];
                        outputGroupSolsInt.append(sol);
                # add profit at dettached output
                intTotal += poolsOutputsAssign[outp][1];
                P = intTotal-initTotal;  # profit of moving out pool/2pools out of old group
                ##############

                ######### Change to new pools for output and attach them to coupled structure
                # go through hierarchy for output outp and find the costs associated with all moves (C)
                movesPerO = []; onlyCostsForO = [];
                for ite, poolComb in enumerate(solPairForO):
                    if ite==3:
                        a=1;
                    #ite = ite+19;poolComb=solPairForO[19];
                    # if not changing pools, change in profit is 0
                    if poolComb[0] == poolsOutputsAssign[outp][0]:
                        movesPerO.append([0, poolComb[0], poolComb[0]]); onlyCostsForO.append(0); continue;
                    if poolComb[0]==[-1]:
                        movesPerO.append([P+poolComb[1] - poolsOutputsAssign[outp][1], poolComb[0], poolComb[0]]);
                        onlyCostsForO.append(0); continue;
                    # if changing pools but remaining at same concentrations, change in profit is just PL.


                    # calculate PL (profit loss of just coupling independently output with new poolComb)
                    #PL = poolComb[1] - poolsOutputsAssign[outp][1]; --- is implicit between new and int Totals
                    PL=0;

                    newTotal =0;
                    poolNew = poolComb[0]; # new pools we (potentially) move to
                    if type(poolNew) is not list: poolNew = [poolNew];
                    poolsOutputsAssignNew = copy.deepcopy(poolsOutputsAssignInt);
                    poolsOutputsAssignNew[outp] = poolComb;
                    seen = [];
                    uniquePairPoolsNew = [item[0] for item in poolsOutputsAssignNew if item[0] != [-1] and
                                len(item[0])==2 and item[0] not in seen and not seen.append(item[0])];
                    coupledPoolListsNew = findPoolGroups(uniquePairPoolsNew);
                    seen = [item for sublist in coupledPoolListsNew for item in sublist];
                    uniqueOnePoolsNew = [item[0][0] for item in poolsOutputsAssignNew if item[0] != [-1] and
                                   len(item[0])==1 and item[0][0] not in seen and not seen.append(item[0][0])];
                    outputGroupSolsNew = [];
                    # add changed coupled pool groups by addition of poolNew
                    for pools in [cp for cp in coupledPoolListsNew if set(cp) & set(poolNew)]:  # for coupled pools sharing any pool with poolNew pair
                        # if poolsChanged and not (set(pools) & set(poolsChanged)) and not (set(poolNew) & set(poolsChanged)):
                        #     C=onlyCostsForAllO_prev[outp][ite];
                        #     onlyCostsForO.append(C);
                        #     total = PL + P + C;
                        #     if abs(total) < 0.0000001: total = 0;
                        #     movesPerO.append([total, poolOld, poolNew]); savedSolves+=1;
                        #     continue;
                        revenue = 0; #revenue2=0;
                        poolsAndSolConcs2 = assocPoolsConcs(poolsOutputsAssignNew);
                        if len(set(x for x in poolsAndSolConcs2 if x[0] in pools)) == len(pools):
                            revenue = sum([x[1] for x in poolsOutputsAssignNew if x[0][0] in pools]);
                        else:
                            outputs = [ix for ix, item in enumerate(poolsOutputsAssignNew) if
                               (item[0][0] in pools or (len(item[0])==2 and item[0][1] in pools))];
                            if len(outputs) == 1 and len(pools) == 2:
                                revenue = PoolO.PoolPairAtOneOutputSol(
                                    outputs[0], pools, allOutputsConcUps, allOutputFlowUps, allOutputProfits,
                                    allInputsConc, allInputCosts, inputsForPools, directsForOutputs, TOL);
                            else:
                                reducedPoolCostMatrix = np.zeros(shape=(len(pools), len(outputs)));
                                for ip, pool in enumerate(pools):
                                    for io, output in enumerate(outputs):
                                        poolAtOut = poolsOutputsAssignNew[output][0];
                                        if pool == poolAtOut or (type(poolAtOut) is list and pool in poolAtOut):
                                            reducedPoolCostMatrix[ip, io] = poolCostMatrix[pool, poolNb + output];
                                probCostMatrix = np.hstack(
                                    [np.array([np.arange(1, len(inputOrigCostMatrix[:, 0]) + len(pools) + 1)]).T,
                                     np.vstack([inputOrigCostMatrix[:, pools], np.zeros(shape=(len(pools), len(pools)))]),
                                     np.vstack([inputOrigCostMatrix[:, poolNb + np.array(outputs)], reducedPoolCostMatrix])]);
                                lines3 = lines2[:];
                                res = formatAndRunGams(lines3, ws, probCostMatrix, pools, outputs,
                                                           inputNb, poolNb, outputNb, allOutputFlowUps, allOutputsConcUps);
                                revenue = res[0];
                        newTotal += revenue;
                        outputGroupSolsNew.append((pools, [revenue]));
                    # add changed one pool groups by addition of poolNew
                    for pool in (set(uniqueOnePoolsNew) - set(uniqueOnePoolsInt)) | (set(uniqueOnePoolsNew) & set(poolNew)):  # for sets of one pool that have been added or that share the pool with the poolNew
                        # if poolsChanged and not (set([pool]) & set(poolsChanged)) and not (set(poolNew) & set(poolsChanged)):
                        #     C=onlyCostsForAllO_prev[outp][ite];
                        #     onlyCostsForO.append(C);
                        #     total = PL + P + C;
                        #     if abs(total) < 0.0000001: total = 0;
                        #     movesPerO.append([total, poolOld, poolNew]); savedSolves+=1;
                        #     continue;
                        outputs = [];
                        for (output, item) in enumerate(poolsOutputsAssignNew):
                            if item[0] == [pool]: outputs.append(output);
                        #### Solve multiple outputs problem in the original group with the output outp removed and independent
                        inputsToPool = inputsForPools[pool];
                        probOrig = cy.I_1_J_unconstr( \
                            # inputs array to pool 1 with its concentrations and costs
                            [allInputsConc[i] for i in inputsToPool], [allInputCosts[i] for i in inputsToPool], \
                            # directs arrays to each output, followed by arrays of all directs concentrations and costs
                            [directsForOutputs[o] for o in outputs], allInputsConc, allInputCosts, \
                            # lower and upper bound concentrations at each output
                            [0] * len(outputs), [allOutputsConcUps[o] for o in outputs], \
                            # upper limit of flow at each output
                            [allOutputFlowUps[o] for o in outputs], \
                            # profit per unit flow sent to each output
                            [allOutputProfits[o] for o in outputs], \
                            # inputs array to pool 1 with its indices
                            inputsToPool, TOL, DEC);
                        # pools=[pool];
                        # reducedPoolCostMatrix = np.zeros(shape=(len(pools), len(outputs)));
                        # for ip, pool in enumerate(pools):
                        #     for io, output in enumerate(outputs):
                        #         poolAtOut = poolsOutputsAssignNew[output][0];
                        #         if pool == poolAtOut or (type(poolAtOut) is list and pool in poolAtOut):
                        #             reducedPoolCostMatrix[ip, io] = poolCostMatrix[pool, poolNb + output];
                        # probCostMatrix = np.hstack(
                        #     [np.array([np.arange(1, len(inputOrigCostMatrix[:, 0]) + len(pools) + 1)]).T,
                        #      np.vstack([inputOrigCostMatrix[:, pools], np.zeros(shape=(len(pools), len(pools)))]),
                        #      np.vstack([inputOrigCostMatrix[:, poolNb + np.array(outputs)], reducedPoolCostMatrix])]);
                        # lines3 = lines2[:];
                        # revenue2 = formatAndRunGams(lines3, ws, probCostMatrix, pools, outputs,
                        #                            inputNb, poolNb, outputNb, allOutputFlowUps, allOutputsConcUps);
                        revenue = probOrig.solve_I_1_J()[0][0];
                        newTotal += revenue;
                        outputGroupSolsNew.append((pool, [revenue]));
                    changed = [cp for cp in coupledPoolListsNew if set(cp) & set(poolNew)];
                    changed = [item for sublist in changed for item in sublist] + \
                              list((set(uniqueOnePoolsNew) - set(uniqueOnePoolsInt)) | (set(uniqueOnePoolsNew) & set(poolNew)));
                    # add unchanged coupled or single pool groups
                    for sol in outputGroupSolsInt:
                        if type(sol[0]) is not tuple:pools = [sol[0]];
                        else: pools = sol[0];
                        if not set(pools) & set(changed) and not set(pools) & set(poolNew):
                            newTotal += sol[1][0];
                            outputGroupSolsNew.append(sol);
                    C = newTotal-intTotal; # cost of coupling poolNew with the relevant structures sharing pools

                    ###### Append to onlyCostsForO buffer, calculate total and add to moves
                    onlyCostsForO.append(C);
                    total = PL + P + C;
                    if abs(total)<0.0000001: total =0;
                    movesPerO.append([total, poolOld, poolNew]);
                    a=1;
                movesForAllO.append(movesPerO);
                bestMovePerO[outp] = max(movesPerO, key=lambda x: x[0]);
                onlyCostsForAllO.append(onlyCostsForO);
            try:
                tradeoffs = [move[0] for move in bestMovePerO];
                outputChanged = tradeoffs.index(max(tradeoffs))
                actualMove = bestMovePerO[outputChanged];
            except TypeError:
                a=1;
            if actualMove[0]<0.000001: #### if no improvement possible then stop
                positiveMove=False;
                filelist = glob.glob( "I:\OneDrive for Business 1\PhD Radu\Ruth-MILP\Python\Pooling Python Code\data\_gams_py*");
                for file in filelist: os.remove(file);
                continue;
            poolOld = actualMove[1]; poolNew = actualMove[2];

            ####################### Perform moving output from one pool to another
            poolsOutputsAssign[outputChanged] = next(poolSol for poolSol in solPairsForAllO[outputChanged] if poolSol[0]==poolNew);

            # calculate updated old output group and new output group
            seen = [];
            uniquePairPools = [item[0] for item in poolsOutputsAssign if item[0] != [-1] and
                               len(item[0]) == 2 and item[0] not in seen and not seen.append(item[0])];
            coupledPoolListsUpd = findPoolGroups(uniquePairPools);

            # find all pools affected directly or indirectly via actual move;
            poolsOld = [list(pools) for pools in coupledPoolLists if set(poolOld)&set(pools)];
            if not poolsOld: poolsOld=[[poolOld]];
            poolsNew = [list(pools) for pools in coupledPoolListsUpd if set(poolNew)&set(pools)];
            if not poolsNew: poolsNew = [[poolNew]];
            poolsChanged = list(itertools.chain.from_iterable(poolsOld+poolsNew));
            poolsChanged = [item for sublist in poolsChanged for item in sublist];

            coupledPoolLists = coupledPoolListsUpd;
            seen = [item for sublist in coupledPoolLists for item in sublist];
            uniqueOnePools = [item[0][0] for item in poolsOutputsAssign if item[0] != [-1] and
                              len(item[0]) == 1 and item[0][0] not in seen and not seen.append(item[0][0])];
            outputGroupSols = [0] * len(uniqueOnePools);
            # solve polynomially all one pool independent groups
            for i, pool in enumerate(uniqueOnePools):
                #### Solve multiple outputs problem in all original groups
                inputsToPool = inputsForPools[pool];
                outputs = [];
                for (output, item) in enumerate(poolsOutputsAssign):
                    if item[0] == [pool]: outputs.append(output);
                # solve one pool, multiple outputs problem
                probOrig = cy.I_1_J_unconstr( \
                    # inputs array to pool 1 with its concentrations and costs
                    [allInputsConc[i] for i in inputsToPool], [allInputCosts[i] for i in inputsToPool], \
                    # directs arrays to each output, followed by arrays of all directs concentrations and costs
                    [directsForOutputs[o] for o in outputs], allInputsConc, allInputCosts, \
                    # lower and upper bound concentrations at each output
                    [0] * len(outputs), [allOutputsConcUps[o] for o in outputs], \
                    # upper limit of flow at each output
                    [allOutputFlowUps[o] for o in outputs], \
                    # profit per unit flow sent to each output
                    [allOutputProfits[o] for o in outputs], \
                    # inputs array to pool 1 with its indices
                    inputsToPool, TOL, DEC);
                outputGroupSols[i] = (pool, probOrig.solve_I_1_J()[0]);
            # solve using gams for any other multi-pool group and their associated outputs
            for pools in coupledPoolLists:
                outputs = [ix for ix, item in enumerate(poolsOutputsAssign) if
                           (item[0][0] in pools or (len(item[0]) == 2 and item[0][1] in pools))]
                reducedPoolCostMatrix = np.zeros(shape=(len(pools), len(outputs)));
                for ip, pool in enumerate(pools):
                    for io, output in enumerate(outputs):
                        poolAtOut = poolsOutputsAssign[output][0];
                        if pool == poolAtOut or (type(poolAtOut) is list and pool in poolAtOut):
                            reducedPoolCostMatrix[ip, io] = poolCostMatrix[pool, poolNb + output];
                probCostMatrix = np.hstack(
                    [np.array([np.arange(1, len(inputOrigCostMatrix[:, 0]) + len(pools) + 1)]).T,
                     np.vstack([inputOrigCostMatrix[:, pools], np.zeros(shape=(len(pools), len(pools)))]),
                     np.vstack([inputOrigCostMatrix[:, poolNb + np.array(outputs)], reducedPoolCostMatrix])]);
                lines3 = lines2[:];
                profitPools = formatAndRunGams(lines3, ws, probCostMatrix, pools, outputs,
                                               inputNb, poolNb, outputNb, allOutputFlowUps, allOutputsConcUps);
                outputGroupSols.append((pools, [profitPools[0]]));
            onlyCostsForAllO_prev = onlyCostsForAllO;

            ##################### Print iteration information and clean temporary files for GAMS
            iterCount = iterCount+1;
            objVal = sum([outputGroupSol[1][0] for outputGroupSol in outputGroupSols]) \
                     + sum([x[1] for x in poolsOutputsAssign if x[0]==[-1]]);
            t1=time.time();
            string = "Iter %d, obj = %.6f, change = +%.6f" % (iterCount, objVal, actualMove[0]) +\
                  "\t\t- pools " + str(poolOld) + "->" + str(poolNew) + "\t\t\tfor output " \
                  + str(outputChanged) + ("\t - time= %f s" % (t1-t0));
            with open(filename, "a") as text_file: text_file.write(string+'\n'); print(string);
            ##### clean temporary gams files
            filelist = glob.glob(
                "I:\OneDrive for Business 1\PhD Radu\Ruth-MILP\Python\Pooling Python Code\data\_gams_py*");
            for file in filelist:
                os.remove(file);
        endTime = time.time();
        string1 = 'Finished iterations for quality # %d - time= %f s' % ((k+1), endTime-startTime);
        string2 = '*RESULTS for quality # %d : ALGO = %.6f BARON lower = %.6f , BARON upper = %.6f' % ((k+1), objVal, lowerB, upperB);
        string3 = '*TIMES   for quality # %d : ALGO = %f s, BARON = %f s' % ((k+1), endTime-startTime, timeBARON);
        with open(filename, "a") as text_file:
            text_file.write('\n' + string1 + '\n');
            text_file.write(string2 + '\n');
            text_file.write(string3 + '\n');
            print(string1);print(string2);print(string3);
    return 0;

        # # group all outputs that receive same concentration at optimality around most common pool they share
        # # no compromise made here
        # optimalSolsForAllO = [];
        # for solPairForO in solPairsForAllO:
        #     optimalSolsForAllO = optimalSolsForAllO + [(elem[0], elem[2]) for elem in solPairForO if
        #                                                elem[1] == solPairForO[0][1] and len(solPairForO) != 1];
        # nbOccurs = Counter(optimalSolsForAllO).most_common();
        # while True:
        #     if nbOccurs and nbOccurs[0][1] > 1:
        #         for i in range(0, len(solPairsForAllO)):
        #             chosen = [elem for elem in solPairsForAllO[i] if
        #                       elem[0] == nbOccurs[0][0][0] and elem[2] == nbOccurs[0][0][1]];
        #             if chosen: solPairsForAllO[i] = chosen;
        #         optimalSolsForAllO = [];
        #         for solPairForO in solPairsForAllO:
        #             optimalSolsForAllO = optimalSolsForAllO + [(elem[0], elem[2]) for elem in solPairForO if
        #                                                        elem[1] == solPairForO[0][1] and len(solPairForO) != 1];
        #         nbOccurs = Counter(optimalSolsForAllO).most_common();
        #     else:
        #         break;
        #
        # # make choice in turn for every output that has an independent pool, separable from the rest
        # optimalSolsForAllO = [];
        # for solPairForO in solPairsForAllO:
        #     optimalSolsForAllO = optimalSolsForAllO + [elem[0] for elem in solPairForO if
        #                                                elem[1] == solPairForO[0][1] and len(solPairForO) != 1];
        # nbOccurs = list(reversed(Counter(optimalSolsForAllO).most_common()));
        # while True:
        #     # make choice in turn for every output that has an independent pool, separable from the rest
        #     if nbOccurs and nbOccurs[0][1] == 1:
        #         for i in range(0, len(solPairsForAllO)):
        #             chosen = [elem for elem in solPairsForAllO[i] if elem[0] == nbOccurs[0][0]];
        #             if chosen: solPairsForAllO[i] = chosen;
        #         optimalSolsForAllO = [];
        #         for solPairForO in solPairsForAllO:
        #             optimalSolsForAllO = optimalSolsForAllO + [elem[0] for elem in solPairForO if
        #                                                        elem[1] == solPairForO[0][1] and len(solPairForO) != 1];
        #         nbOccurs = list(reversed(Counter(optimalSolsForAllO).most_common()));
        #     # from among those outputs with pools undecided, decide them in turn starting with the one where going to second-best option costs most
        #     elif nbOccurs:
        #         costsSecondBest = [];
        #         for i in range(0, len(solPairsForAllO)):
        #             if len(solPairsForAllO[i]) > 1:
        #                 costsSecondBest.append([i, solPairsForAllO[i][0][1] - solPairsForAllO[i][1][1]]);
        #             else:
        #                 costsSecondBest.append([i, 0]);
        #         # take output with most to lose from choosing next best pools for it and see if that or combo of first choice is best
        #         [indexOutput, costSecondChoice] = max(costsSecondBest, key=lambda x: x[1]);
        #         poolFirstChoice = solPairsForAllO[indexOutput][0][0];
        #         outputsCoupled = [index for (index, solPairsForO) in enumerate(solPairsForAllO)
        #                           if len(solPairsForO) == 1 and solPairsForO[0][0] == 3];
        #         inputsToPool = inputsForPools[poolFirstChoice];
        #         # construct coupling solution and compare with second choice solution in terms of loss
        #         prob = cy.I_1_J_unconstr( \
        #             # inputs array to pool 1 with its concentrations and costs
        #             [allInputsConc[i] for i in inputsToPool], [allInputCosts[i] for i in inputsToPool], \
        #             # directs arrays to each output, followed by arrays of all directs concentrations and costs
        #             [directsForOutputs[o] for o in outputsCoupled], allInputsConc, allInputCosts, \
        #             # lower and upper bound concetrations at each output
        #             [0] * len(outputsCoupled), [allOutputsConcUps[o] for o in outputsCoupled], \
        #             # upper limit of flow at each output
        #             [allOutputFlowUps[o] for o in outputsCoupled], \
        #             # profit per unit flow sent to each output
        #             [allOutputProfits[o] for o in outputsCoupled], \
        #             # inputs array to pool 1 with its indices
        #             inputsToPool, TOL, DEC);
        #         ans = prob.solve_I_1_J;
        #
        #         a = 1;
        #     else:
        #         break;  # if all outputs have chosen pool


        # allOutputPoolsComb = list(itertools.product(*pairsForAllO));
        # adjMatrixModif  = np.zeros((inputNb+poolNb,  poolNb+outputNb+1));
        # for idx,lineText in enumerate(lines2):

# if "a(i,j)" in line:
#                inputAdjRows = lines[idx+2         : idx+2+inputNb];
#                poolAdjRows  = lines[idx+2+inputNb : idx+2+inputNb+poolNb];
#                for i,l in enumerate(poolAdjRows):
#                    lines[idx+2+inputNb+i]=" ";
#                    adjArrayForPool  = np.fromstring(l, dtype=float, sep=' ')[1:];
#
#
#
#
#    if len(sys.argv) > 1:
#        ws = GamsWorkspace(system_directory = sys.argv[1])
#    else:
#        ws = GamsWorkspace()
#    ws = GamsWorkspace(debug=DebugLevel.KeepFiles)
#    #file = open(os.path.join(ws.working_directory, "tdata.gms"), "a")
#    #file.write(get_data_text())
#    #file.close()
#    #file = open(os.path.join(ws.working_directory, "Haverly1.gms"), "r")
#    data = '';
#    with open(os.path.join(os.path.curdir,problemName), "r") as dataFile:
#        data = dataFile.read();
#    #file = open(os.path.join(ws.working_directory, "Haverly1.gms"), "r")
#    t2 = ws.add_job_from_string(data);
#    opt = ws.add_options()
#    #opt.defines["incname"] = "pqmodel"
#    t2.run(opt)
#    for rec in t2.out_db["f"]:
#        print("f(" + rec.key(0) + "," + rec.key(1) + "): level=" + str(rec.level))

    # prob2 = cy.I_1_J_unconstr( \
    #     # inputs array to pool 1 with concentrations and costs
    #     [3, 1, 4, 2, 1.5], [6, 13, 10, 9, 11], \
    #     # directs arrays to each output, followed by arrays of all directs concentrations and costs
    #     [[0, 1, 2], [1, 2]], [1.5, 3, 5], [10, 11, 9.5], \
    #     # lower and upper bound concetrations at each output
    #     [2, 1.5], [3.5, 3], \
    #     # upper limit of flow at each output
    #     [200, 250], \
    #     # profit per unit flow sent to each output
    #     [15, 14]);
    #
    # prob1 = cy.I_1_J_unconstr( \
    #     # inputs array to pool 1 with concentrations and costs
    #     [3, 1], [6, 13], \
    #     # directs arrays to each output, followed by arrays of all directs concentrations and costs
    #     [[0], [1]], [1.5, 2.5], [12, 3], \
    #     # lower and upper bound concetrations at each output
    #     [1, 1], [2, 2], \
    #     # upper limit of flow at each output
    #     [200, 150], \
    #     # profit per unit flow sent to each output
    #     [15, 15]);

main();
# ans = prob1.solve_I_1_1(1); # solve for output 1
# ans = solve_I_1_J(prob2,tol,dec); # solve 2 output problem
# print(ans)



# prob1 = I_1_1_unconstr(tol, \
#    # inputs array to pool 1 with concentrations and costs
#    [0,1,2,3,4], [3,1,4,2,1.5], [6,13,10,9,11], \
#    # directs arrays to each output, followed by arrays of all directs concentrations and costs
#    [0,1,2],     [1.5, 3, 5],   [10, 11, 9.5], \
#    # lower and upper bound concetrations at each output
#    [2], [3.5], \
#    # upper limit of flow at each output
#    [200], \
#    # profit per unit flow sent to each output
#    [15]);
# ans = solve_I_1_1(prob1,tol,dec);

# *****************    CORRECTIONS MADE FROM WRITEUP    ************************
#    - formula for finding breakpoint roots between two mixed sets was slightly wrong, corrected it in the code
#    - extra condition needed when choosing between P_L/P_U as a mixed concentration for mixed triples, included in code
#    - needed directs vs either inputs or mixed breakpoints for the multiple outputs case,
#        need to include formuals/derivations in writeup

# ****************     THINGS TO STILL CONSIDER        *************************
#    - code refactoring to not pass input object around and more code reuse
#    - sometimes an identical solution can be found with two flow configurations (atm I am just choosing a random one)
#        so need to consider alternatives which are important in the constrained case
#    - discontinuities can arise at a breakpoint (not just differentiability breaks down)
#        so need to consider active sets on both sides (for 1 output problems)
#    - extension to find exact constrained flow solution to I_1_J problem:
#        *need per unit flow comparison of profit/cost
#        *need full domination ranking of all input pairs/mixed triples
#        *need updating mechanism of bounds after every sweep

