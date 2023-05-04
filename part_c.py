import gurobipy as gp
import time

s = "AAGUUUUAGGAGCAGCCUUAGUGUGAACAGCAAUGCCAUAAUUGAGUCACAAGUUGCCAAAGCGUGACAU"
startTime = time.time()

# model
model = gp.Model('RNA_Pairing')

# Define the decision variables
n = len(s)
X = model.addVars(n, n, vtype=gp.GRB.BINARY, name='X')

# objective function
model.setObjective(gp.quicksum(X[i, j] * (-1.33 if s[i] == 'A' or s[i] == 'U' else -1.45) for i in range(n) for j in range(i+1, n)), gp.GRB.MINIMIZE)

# constraints
for i in range(n):
    model.addConstr(gp.quicksum(X[i, j] + X[j, i] for j in range(i+1, n)) <= 1)
    model.addConstr(gp.quicksum(X[i, j] for j in range(i+1, n)) <= 1)
for j in range(n):
    model.addConstr(gp.quicksum(X[i, j] for i in range(0, j)) <= 1)

for i in range(n):
    for j in range(i, n):
        if j >= i+7:
            model.addConstr(X[i, j] <= 1)
        else:
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
                print(f"Pair ({i},{j}), {s[i]} : {s[j]}")
                c += 1
    print(f"Optimal solution found: {c} pairs")
    print(f"energy is {model.objVal}")
    print(f"elapsed time is: {time.time() - startTime}")
else:
    print("No feasible solution found.")
