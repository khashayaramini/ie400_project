import gurobipy as gp
import pandas as pd
import time

dataFile = 'data.xlsx'
df = pd.read_excel(dataFile)
vals = df.iloc[:, 1].values

print(vals[5])
s = vals[5]

keys = [AU, CG, GC, GU, UG, UA]
energies = [[-1.1, -2.1, -2.2, -1.4, -0.9, -0.6],[-2.1, -2.4, -3.3, -2.1, -2.1, -1.4],[-2.2, -3.3, -3.4, -2.5, -2.4, -1.5],[-1.4, -2.1, -2.5, -1.3, -1.3, -0.5],[-0.9, -2.1, -2.4, -1.3, -1.3, -1.0],[-0.6, -1.4, -1.5, -0.5, -1.0, -0.3]]

startTime = time.time()

# model
model = gp.Model('RNA_Pairing')

# Define the decision variables
n = len(s)
X = model.addVars(n, n, vtype=gp.GRB.BINARY, name='X')
Y = model.addVars(n, n, vtype=gp.GRB.BINARY, name='Y')

# objective function
model.setObjective(gp.quicksum(X[i, j] for i in range(n) for j in range(i+1, n)), gp.GRB.MAXIMIZE)
model.setObjective(gp.quicksum(X[i, j] for i in range(n) for j in range(i+1, n)), gp.GRB.MAXIMIZE)

# constraints

for i in range(n):
	model.addConstr((gp.quicksum(Y[i, j] + Y[i+1, j-1] ) == 2)
	
for i in range(n):
    model.addConstr(gp.quicksum(X[i, j] + X[j, i] for j in range(i+1, n)) <= 1)
    model.addConstr(gp.quicksum(X[i, j] for j in range(i+1, n)) <= 1)
for j in range(n):
    model.addConstr(gp.quicksum(X[i, j] for i in range(0, j)) <= 1)

for i in range(n):
    for j in range(i, n):
        if j >= i+4:
            model.addConstr(X[i, j] <= 1)
        if j == i+1 or j == i+2 or j == i+3:
            model.addConstr(X[i, j] == 0)
        if X[i, j]:
            if s[i] == 'A' and s[j] == 'U' or s[i] == 'U' and s[j] == 'A':
                model.addConstr(X[i, j] <= 1)
            elif s[i] == 'G' and s[j] == 'C' or s[i] == 'C' and s[j] == 'G':
                model.addConstr(X[i, j] <= 1)
            else:
                model.addConstr(X[i, j] == 0)
                model.addConstr(X[j, i] == 0)
	
for i in range(n):
    for j in range(i+1, n):
        for k in range(i+1, j):
            for l in range(j+1, n):
                model.addConstr(X[i, j]*X[k, l] + X[i, k]*X[j, l] + X[i, l]*X[j, k] == 0)

model.Params.OutputFlag = 0  # suppress Gurobi output

	
# Optimize
model.optimize()

# printing result
c = 0
if model.status == gp.GRB.OPTIMAL:
    for i in range(n):
        for j in range(i+1, n):
            if round(X[i, j].x) == 1:
                print(f"Pair ({i+1},{j+1}), {s[i]} : {s[j]}")
                c += 1
    print(f"Optimal solution found: {c} pairs")
    print(f"elapsed time is: {time.time() - startTime}")
else:
    print("No feasible solution found.")
