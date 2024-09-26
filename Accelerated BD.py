from pyomo.environ import *
import numpy as np
import pandas as pd
from datetime import datetime

def delete_component(Model, comp_name):

    list_del = [vr for vr in vars(Model)
                if comp_name == vr
                or vr.startswith(comp_name + '_index')
                or vr.startswith(comp_name + '_domain')]

    for k in list_del:
        Model.del_component(k)

#%% #----------Excel import----------#
address=r'C:\Users\Amirreza\Desktop\Amirreza\Omega\Small Data\S14.xlsx'
node_data = pd.read_excel(address, sheet_name='Node')
Node = pd.DataFrame(node_data)
edge_data = pd.read_excel(address, sheet_name='Edge')
Edge = pd.DataFrame(edge_data)
refugee_data = pd.read_excel(address, sheet_name='Refugee')
Refugee = pd.DataFrame(refugee_data)
MF_data = pd.read_excel(address, sheet_name='Facility')
Facility = pd.DataFrame(MF_data)
distance_data = pd.read_excel(address, sheet_name='Distance')
Distance = pd.DataFrame(distance_data)
travel_data = pd.read_excel(address, sheet_name='Travel')
Travel = pd.DataFrame(travel_data)

#%% #----------Sets & Parameters----------#
R = list(Refugee['R'])
p = {}
e = {}
d = {}
for r in R:
    p[r] = Refugee['p'][r-1]
    e[r] = Refugee['e'][r-1]
    d[r] = Refugee['d'][r-1]
M = list(Facility['M'])
M2 = M[:-1]
b = {}
f = {}
o = {}
for m in M:
    b[m] = Facility['b'][m-1]
    f[m] = Facility['f'][m-1]
    o[m] = Facility['o'][m-1]
VR = list(Node['Vr'])
VM = list(Node['Vm'])
VRM = list(Node['Vrm'])
P = list(Edge['P'])
l = {}
for path in P:
    l[path] = Edge['l'][path-1]
K = max(l.values())
n = {}
for path in P:
    for k in RangeSet(K):
        if Edge[k][path-1] > 0 : 
            n[path,k] = Edge[k][path-1]
tau = int(Refugee['tau'][0])
LR = {}
Tmax = 0
for r in R:
    Lrlist = []
    for q in RangeSet(0 , l[p[r]] - 1):
        Lrlist.append((n[p[r],q+1] , e[r]+q))
        if q == l[p[r]] - 1:
            if e[r]+q > Tmax:
                Tmax = e[r]+q
    LR[r] = Lrlist
T = RangeSet(Tmax)
L = []
for r in R:
    for j in LR[r]:
        if j not in L:
            L.append(j)
J = {}
c = {}
for i in VM:
    acc = []
    for j in VM:
        if Distance[j][i-1] <= 200:
            acc.append(j)
            c[i,j] = Travel[j][i-1]
    J[i] = acc

#%% #----------Model----------#
A_IDX = []
for r in R:
    for i,t in LR[r]:
        A_IDX.append((r,i,t))
Y_IDX = []
for m in M:
    for i,t in L:
        Y_IDX.append((m,i,t))
X_IDX = []
for m in M:
    for i in VM:
        for j in J[i]:
            for t in T:
                X_IDX.append((m,i,j,t))
Z_IDX = []
for m in M:
    Z_IDX.append(m)

#%% #----------CM----------#
# Storing matrices
CM_Alist = {}
for idx in A_IDX:
    CM_Alist[idx] = 0
CM_Ylist = {}
for idx in Y_IDX:
    CM_Ylist[idx] = 0
CM_Xlist = {}
for idx in X_IDX:
    CM_Xlist[idx] = 0
CM_Zlist = {}
for idx in Z_IDX:
    CM_Zlist[idx] = 0
solver = SolverFactory('gurobi')
sorted_P = list(dict(sorted(l.items(), key=lambda x:x[1])).keys()) # Sort paths in a non-increasing order of length
R_bar = []  # Set of RGs included in partial network N
R_old = []  # RGs inserted in the previous iterations.
iterr = 0
runTlist = []
print('#-------------Welcome to CM-------------#')
CM_X_IDX = []
for m in M:
    for i in VM:
        for j in J[i]:
            for t in T:
                CM_X_IDX.append((m,i,j,t))
CM_Z_IDX = []
for m in M:
    CM_Z_IDX.append(m)

for p_index in sorted_P:  # Main loop
    iterr += 1
    R_new = [r for r in R if p[r] == p_index]
    for r in R_new:
        R_bar.append(r)
    CM_L = []
    for r in R_bar:
        for j in LR[r]:
            if j not in CM_L:
                CM_L.append(j)
    # CM indices:
    CM_A_IDX = []
    for r in R_bar:
        for i,t in LR[r]:
            CM_A_IDX.append((r,i,t))
    CM_Y_IDX = []
    for m in M:
        for i,t in CM_L:
            CM_Y_IDX.append((m,i,t))
    
    # Model
    CM = ConcreteModel()
    # Variables
    CM.A = Var(CM_A_IDX, within = NonNegativeReals)
    CM.Y = Var(CM_Y_IDX, within = Binary)
    CM.X = Var(CM_X_IDX, within = Binary)
    CM.Z = Var(CM_Z_IDX, within = Binary)

    CM.obj = Objective(expr = sum(f[m]*CM.Z[m] for m in CM_Z_IDX) + 
                   sum(o[m]*CM.Y[m,i,t] for m in M for i,t in CM_L) +
                   sum(c[i,j]*CM.X[m,i,j,t] for m in M for i in VM for j in J[i] for t in T))

    CM.limit = ConstraintList()
    for r in R_bar:
        for q in RangeSet(0,l[p[r]]-tau):
            CM.limit.add(sum(CM.A[r,n[p[r],q+t2+1],e[r]+q+t2] for t2 in range(tau))>=1)

    for i,t in CM_L:
        CM.limit.add(sum(b[m] * CM.Y[m,i,t] for m in M)>=sum(d[r]*CM.A[r,i,t] for r in R_bar if (i,t) in LR[r]))

    for m in M:
        for i,t in CM_L:
            CM.limit.add(CM.Y[m,i,t] <= sum(CM.X[m,i,j,t] for j in J[i]))

    for m in M:
        for i in VM:
            for t in RangeSet(Tmax-1):
                CM.limit.add(sum(CM.X[m,j,i,t] for j in VM if i in J[j]) <= sum(CM.X[m,i,j,t+1] for j in J[i]))

    for m in M:
        for t in T:
            CM.limit.add(sum(CM.X[m,i,j,t] for i in VM for j in J[i]) <= CM.Z[m])

    for m in M2:
        CM.limit.add(CM.Z[m]>=CM.Z[m+1])
    
    # Use of stored data
    if iterr > 1:
        for idx in CM_A_IDX:
            if idx[0] in R_old:
                CM.limit.add(CM.A[idx] == CM_Alist[idx])
        for idx in CM_Y_IDX:
            if CM_Ylist[idx] > 0.5:
                CM.limit.add(CM.Y[idx] == 1)
    start_time = datetime.now()
    solver.solve(CM) 
    end_time = datetime.now()
    runT = end_time - start_time
    runTlist.append(runT.seconds)
    print('-----Iteration', iterr,'-----')
    print('Path inserted:', p_index)
    print('RGs included:', R_new)
    print('Objective:', int(round(value(CM.obj))))
    print('Partial run time:', runT.seconds)
    for r in R_new:
        R_old.append(r)
    # Storing the data
    for idx in CM_A_IDX:
        CM_Alist[idx] = value(CM.A[idx])
    for idx in CM_Y_IDX:
        CM_Ylist[idx] = value(CM.Y[idx])

for idx in CM_X_IDX:
    CM_Xlist[idx] = value(CM.X[idx])
for idx in CM_Z_IDX:
    CM_Zlist[idx] = value(CM.Z[idx])

#print(address[-8:-5])
print('#-------------CM overall results-------------#')
print('Total CM objective:', int(round(value(CM.obj))))
print('Total MFs utilized:', int(value(sum(CM.Z[m] for m in CM_Z_IDX))))
print('Total CM run-time:', sum(runTlist))

#%% #----------Model----------#
A_IDX = []
for r in R:
    for i,t in LR[r]:
        A_IDX.append((r,i,t))
Y_IDX = []
for m in M:
    for i,t in L:
        Y_IDX.append((m,i,t))
X_IDX = []
for m in M:
    for i in VM:
        for j in J[i]:
            for t in T:
                X_IDX.append((m,i,j,t))
Z_IDX = []
for m in M:
    Z_IDX.append(m)

#%% #-------------------------DSP-------------------------#
DSP = ConcreteModel()
# PArams
Ybar = {}
for i in Y_IDX:
    Ybar[i] = 0

u_IDX = []
for r in R:
    for q in RangeSet(0,l[p[r]] - tau):
        u_IDX.append((r,q))
v_IDX = []
for i in L:
    v_IDX.append(i)
#%%
# Variables
DSP.u = Var(u_IDX, within = NonNegativeReals)
DSP.v = Var(v_IDX, within = NonNegativeReals)
DSP.ybar = Param(Y_IDX, initialize = 0, mutable = True)
B = 10000
# Objective
DSP.obj = Objective(expr = sum(DSP.u[r,q] for r, q in u_IDX) - sum(b[m] * DSP.ybar[(m,i,t)]*DSP.v[i,t] for m in M for i, t in v_IDX)
                    , sense=maximize)

# Constraints
DSP.limit = ConstraintList()
for r in R:
    for k in RangeSet(l[p[r]]):
        a1 = max(0 , k - tau)
        a2 = min(k - 1 , l[p[r]] - tau)
        DSP.limit.add(sum(DSP.u[r,q] for q in RangeSet(a1,a2)) - d[r] * DSP.v[n[(p[r],k)],e[r]+k-1] <= 0)

DSP.limit.add(sum(DSP.u[r,q] for r, q in u_IDX) - sum(b[m] * DSP.ybar[(m,i,t)]*DSP.v[i,t] for m in M for i, t in v_IDX) <= B)
DSP.temp = ConstraintList()
#%% 
RMP = ConcreteModel()
# Variables
RMP.Y = Var(Y_IDX, initialize = 0, within = Binary)
RMP.X = Var(X_IDX, initialize = 0, within = Binary)
RMP.Z = Var(Z_IDX, initialize = 0, within = Binary)
ubar = {}
for idx in u_IDX:
    ubar[idx] = 0
vbar = {}
for idx in v_IDX:
    vbar[idx] = 0

RMP.obj = Objective(expr = sum(f[m]*RMP.Z[m] for m in Z_IDX) + 
                   sum(o[m]*RMP.Y[m,i,t] for m in M for i,t in L) +
                   sum(c[i,j]*RMP.X[m,i,j,t] for m in M for i in VM for j in J[i] for t in T))

RMP.limit = ConstraintList()
for m in M:
    for i,t in L:
            RMP.limit.add(RMP.Y[m,i,t] <= sum(RMP.X[m,i,j,t] for j in J[i]))

for m in M:
    for i in VM:
        for t in RangeSet(Tmax-1):
            RMP.limit.add(sum(RMP.X[m,j,i,t] for j in VM if i in J[j]) <= sum(RMP.X[m,i,j,t+1] for j in J[i]))

for m in M:
    for t in T:
        RMP.limit.add(sum(RMP.X[m,i,j,t] for i in VM for j in J[i]) <= RMP.Z[m])

for m in M2:
    RMP.limit.add(RMP.Z[m]>=RMP.Z[m+1])
RMP.limit.add(sum(RMP.Z[m] for m in Z_IDX) >= ceil(sum(d.values()) / ((sum(b.values()) / len(b)) * tau)))

RMP.cut = ConstraintList()
for r, q in u_IDX:
    RMP.cut.add(0 >= d[r] - sum(b[m] * RMP.Y[m,LR[r][q+k-1][0],LR[r][q+k-1][1]] for m in M for k in RangeSet(tau)))

#%% LRRMP
LRRMP = ConcreteModel()
# Variables
LRRMP.Y = Var(Y_IDX, within = NonNegativeReals, bounds = (0,1))
LRRMP.X = Var(X_IDX, within = NonNegativeReals, bounds = (0,1))
LRRMP.Z = Var(Z_IDX, within = NonNegativeReals, bounds = (0,1))
LRRMP.ubar = Param(u_IDX, initialize = 0)

LRRMP.obj = Objective(expr = sum(f[m]*LRRMP.Z[m] for m in Z_IDX) + 
                   sum(o[m]*LRRMP.Y[m,i,t] for m in M for i,t in L) +
                   sum(c[i,j]*LRRMP.X[m,i,j,t] for m in M for i in VM for j in J[i] for t in T))

LRRMP.limit = ConstraintList()
for m in M:
    for i,t in L:
            LRRMP.limit.add(LRRMP.Y[m,i,t] <= sum(LRRMP.X[m,i,j,t] for j in J[i]))

for m in M:
    for i in VM:
        for t in RangeSet(Tmax-1):
            LRRMP.limit.add(sum(LRRMP.X[m,j,i,t] for j in VM if i in J[j]) <= sum(LRRMP.X[m,i,j,t+1] for j in J[i]))

for m in M:
    for t in T:
        LRRMP.limit.add(sum(LRRMP.X[m,i,j,t] for i in VM for j in J[i]) <= LRRMP.Z[m])

for m in M2:
    LRRMP.limit.add(LRRMP.Z[m]>=LRRMP.Z[m+1])
LRRMP.limit.add(sum(LRRMP.Z[m] for m in Z_IDX) >= ceil(sum(d.values()) / ((sum(b.values()) / len(b)) * tau)))

LRRMP.cut = ConstraintList()
for r, q in u_IDX:
    LRRMP.cut.add(0 >= d[r] - sum(b[m] * LRRMP.Y[m,LR[r][q+k-1][0],LR[r][q+k-1][1]] for m in M for k in RangeSet(tau)))

#%%
switch = False
maxIt = 30
iterr = 0

while iterr<maxIt:      # Iterations
    print("-----Iteration",iterr,"-----")
    if iterr==1:
        start_time = datetime.now()
    if switch == False:
        solver.solve(LRRMP, timelimit = 21600)
        LB = value(LRRMP.obj)
        for idx in Y_IDX:
            DSP.ybar[idx]=value(LRRMP.Y[idx])
    else:
        print('#-----The feasible solution of CM is fed to the RMP-----#')
        for idx in Y_IDX:
            RMP.Y[idx] = round(CM_Ylist[idx])   
        for idx in X_IDX:
            RMP.X[idx] = round(CM_Xlist[idx])
        for idx in Z_IDX:
            RMP.Z[idx] = round(CM_Zlist[idx])
        solver.solve(RMP, timelimit = 21600, warmstart = True, options_string="mipgap=0.0001")
        LB = value(RMP.obj)
        for idx in Y_IDX:
            DSP.ybar[idx]=value(RMP.Y[idx])
        #print("RMP obj = ", value(RMP.obj))
    
    print("LB = ", LB)
    
    #print('Eta =', value(RMP.eta))
    
    #for idx in Y_IDX:
    #    DSP.ybar[idx]=value(RMP.Y[idx])

    solver.solve(DSP)
    if switch == False:
        UB = value(sum(f[m]*LRRMP.Z[m] for m in Z_IDX) + 
                   sum(o[m]*LRRMP.Y[m,i,t] for m in M for i,t in L) +
                   sum(c[i,j]*LRRMP.X[m,i,j,t] for m in M for i in VM for j in J[i] for t in T)) + value(DSP.obj)
    else:
        UB = value(sum(f[m]*RMP.Z[m] for m in Z_IDX) + 
                   sum(o[m]*RMP.Y[m,i,t] for m in M for i,t in L) +
                   sum(c[i,j]*RMP.X[m,i,j,t] for m in M for i in VM for j in J[i] for t in T)) + value(DSP.obj)

    print("UB = ", UB) 
    if abs(UB-LB)<=1:
        if switch == False:
            switch = True
        else:
            break

    while True:
        if value(DSP.obj)<=0.01:
            break
        
        for r, q in u_IDX:
            ubar[(r,q)]=value(DSP.u[(r,q)])
        for i, t in v_IDX:
            vbar[(i,t)]=value(DSP.v[(i,t)])
        # Cut
        RMP.cut.add(0 >= sum(ubar[(r,q)] for r, q in u_IDX) - sum(b[m] * vbar[(i,t)] * RMP.Y[m,i,t] for m, i, t in Y_IDX))

        if switch == False:
            LRRMP.cut.add(0 >= sum(ubar[(r,q)] for r, q in u_IDX) - sum(b[m] * vbar[(i,t)] * LRRMP.Y[m,i,t] for m, i, t in Y_IDX))
        
        for r, q in u_IDX:
            if value(DSP.u[(r,q)]) > 0.01:
                DSP.temp.add(DSP.u[(r,q)] == 0)
        solver.solve(DSP)
        
    delete_component(DSP, 'temp')
    DSP.temp = ConstraintList()
    iterr+=1

final_time = datetime.now()
time = final_time - start_time
#%%
print('instance:', address[-8:-5])
print('#-------------CM overall results-------------#')
print('Total CM objective:', int(round(value(CM.obj))))
print('Total MFs utilized:', int(value(sum(CM.Z[m] for m in CM_Z_IDX))))
print('Total CM run-time:', sum(runTlist))

print('       ')
print('#-----------Acc BD Results-----------#')
print('ACC BD RMP optimal obj:', int(value(RMP.obj)))
print('Total MFs utilized:', int(value(sum(RMP.Z[m] for m in Z_IDX))))
print('ACC BD time:', time.seconds, 'seconds')
print('ACC BD iter:', iterr + 1)