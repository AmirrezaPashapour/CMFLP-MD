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
address=r'C:\Users\Amirreza\Desktop\Amirreza\Omega\Small Data\S13.xlsx'
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

#%% #----------Domains----------#
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

print('#-------------CM overall results-------------#')
print('Total CM objective:', int(round(value(CM.obj))))
print('Total MFs utilized:', int(value(sum(CM.Z[m] for m in CM_Z_IDX))))
print('Total CM run-time:', sum(runTlist))

#%% #----------MILP----------#
print('#-------------Welcome to MILP-------------#')
Model = ConcreteModel()
# Variables
Model.A = Var(A_IDX, within = NonNegativeReals)
Model.Y = Var(Y_IDX, within = Binary)
Model.X = Var(X_IDX, within = Binary)
Model.Z = Var(Z_IDX, within = Binary)

Model.obj = Objective(expr = sum(f[m]*Model.Z[m] for m in Z_IDX) + 
                   sum(o[m]*Model.Y[m,i,t] for m in M for i,t in L) +
                   sum(c[i,j]*Model.X[m,i,j,t] for m in M for i in VM for j in J[i] for t in T))

Model.limit = ConstraintList()
for r in R:
    for q in RangeSet(0,l[p[r]]-tau):
        Model.limit.add(sum(Model.A[r,n[p[r],q+t2+1],e[r]+q+t2] for t2 in range(tau))>=1)

for i,t in L:
        Model.limit.add(sum(b[m] * Model.Y[m,i,t] for m in M)>=sum(d[r]*Model.A[r,i,t] for r in R if (i,t) in LR[r]))

for m in M:
    for i,t in L:
            Model.limit.add(Model.Y[m,i,t] <= sum(Model.X[m,i,j,t] for j in J[i]))

for m in M:
    for i in VM:
        for t in RangeSet(Tmax-1):
            Model.limit.add(sum(Model.X[m,j,i,t] for j in VM if i in J[j]) <= sum(Model.X[m,i,j,t+1] for j in J[i]))

for m in M:
    for t in T:
        Model.limit.add(sum(Model.X[m,i,j,t] for i in VM for j in J[i]) <= Model.Z[m])

for m in M2:
    Model.limit.add(Model.Z[m]>=Model.Z[m+1])

to_warm_start = True
if to_warm_start == True:
    print('#-----The feasible solution of CM is fed to the MILP-----#')

    for idx in A_IDX:
        Model.A[idx] = CM_Alist[idx]
    for idx in Y_IDX:
        Model.Y[idx] = round(CM_Ylist[idx])   
    for idx in X_IDX:
        Model.X[idx] = round(CM_Xlist[idx])
    for idx in Z_IDX:
        Model.Z[idx] = round(CM_Zlist[idx])
#-----
start_time = datetime.now()
solver.solve(Model, tee=True, warmstart = to_warm_start, options_string="mipgap=0.0001")
end_time = datetime.now()
runT = end_time - start_time

print('#-------------CM overall results-------------#')
print('Total CM objective:', int(round(value(CM.obj))))
print('Total MFs utilized:', int(value(sum(CM.Z[m] for m in CM_Z_IDX))))
print('Total CM run-time:', sum(runTlist))

print('#-------------MILP overall results-------------#')
print("Optimal MILP objective = ", round(value(Model.obj)))
print("Total MFs utilized:",int(value(sum(Model.Z[m] for m in Z_IDX))))
print('Total MILP run-time:', runT.seconds)