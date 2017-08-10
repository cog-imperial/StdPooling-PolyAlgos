import math
from numpy import isreal
import numpy as np
from scipy import poly1d, roots
from sympy import symbols, Poly
# from gams import GamsWorkspace
# import os
# import sys
# from collections import Counter
#import ctypes as ct
#import pyximport; pyximport.install()

# try:
#     from line_profiler import LineProfiler
#
#
#     def do_profile(follow=[]):
#         def inner(func):
#             def profiled_func(*args, **kwargs):
#                 try:
#                     profiler = LineProfiler()
#                     profiler.add_function(func)
#                     for f in follow:
#                         profiler.add_function(f)
#                     profiler.enable_by_count()
#                     return func(*args, **kwargs)
#                 finally:
#                     profiler.print_stats()
#
#             return profiled_func
#
#         return inner
#
# except ImportError:
#     def do_profile(follow=[]):
#         "Helpful if you accidentally leave in production!"
#
#         def inner(func):
#             def nothing(*args, **kwargs):
#                 return func(*args, **kwargs)
#
#             return nothing
#
#         return inner


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
        return [x for x in XZ_sol if x[0] == maxProfit];#XZ_sol;#

    def solve_I_1_J(self):
        tol = self.tol;
        dec = self.dec;
        Cz = self.Cz;
        activeForAllO = [];
        bpsCommon = [];
        for o in range(0, len(self.TZ)):
            # [Z_sol, Z_Cost] = self.find_Z_sol(o,1);
            XZ_active = self.find_XZ_active(o, 1);  # get all active sets for output o together with their p-derivatives
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
                        active_set[2]] / tol]);  # append derivative factor plus Cq - still need to divide by (p-Cq)^2
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
                    zIndices = [j]; zFlows = [zj / tol];
                elif zj == 0:
                    zIndices = [i]; zFlows = [zi / tol];
                else:
                    zIndices = [i, j]; zFlows = [zi / tol, zj / tol];
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
                    XZ_active.append(X); continue;  # input pair dominates entire interval
                else:
                    R = Q;
            else:
                if Iu > P_U:
                    Q = [z for z in TZ if Cz[z] < P_U];  # mixed triple 2
                else:
                    Q = [z for z in TZ if Cz[z] > P_L];  # mixed triple 3
                if not Q:
                    XZ_active.append(['infeasible', Il, Iu]); continue;  # input/mixed sets infeasible
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
                            S = [];break;
                    else:
                        if b < S[1]:
                            S[0] = b;
                        else:
                            S = [];break;
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
                        # potential breakpoint between mixed triple and direct-only set
                        pNum = (gq * Pq - gZ * Cq) * (Ci - Cj) + (gj * Ci - gi * Cj) * (Cq - Pq);
                        pDenom = (gq - gZ) * (Ci - Cj) - (gi - gj) * (Cq - Pq);
                        p = round(pNum / pDenom, dec);
                        # cost for mixed triple at bp (lower end of interval)
                        gXZAtbp = (gq * (Pq - bp) + gXZAtbp * (Cq - Pq)) / (Cq - bp);
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
                                       d - (Pq - Cq) * (gi * (bp - Cj) - gj * (bp - Ci)) / ((Ci - Cj) * (bp - Cq))) > 0;
                # check profitability conditions for either input or mixed pair
                if bp * pDenom < pNum < bp2 * pDenom and pDenom != 0:  # if there is a valid profitability breakpoint
                    p = pNum / pDenom;
                    if isProfitableAtBp:
                        XZ_active4.append([active_set, bp, p, deriv]);  XZ_active4.append([['infeasible'], p, bp2, []]);
                    else:
                        XZ_active4.append([['infeasible'], bp, p, []]); XZ_active4.append([active_set, p, bp2, deriv]);
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
                p0=P[0]; p1=P[1];
                if p0 != i and p1 != i and (p-Cx[p0])*(p-Cx[p1])<0: # ((Cx[P[0]] <= p <= Cx[P[1]]) or (Cx[P[0]] >= p >= Cx[P[1]])):
                    # gAgreg = (gx[P[0]]*(p-Cx[P[1]]) + gx[P[1]]*(Cx[P[0]]-p)) / (Cx[P[0]]-Cx[P[1]]);
                    gAgregNom = gx[p0] * (p - Cx[p1]) + gx[p1] * (Cx[p0] - p);
                    gAgregDenom = Cx[p0] - Cx[p1];
                    if (gAgregNom < gx[i] * gAgregDenom) * (gAgregDenom > 0)or \
                        (gAgregNom > gx[i] * gAgregDenom) * (gAgregDenom < 0) : isBp = False; continue;
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
            if len(P0) == 1: return [[Xpair + P0, pl, pu]];  # if no breakpoint occurs add only dominating mixed pair
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
# @do_profile(follow=[])
