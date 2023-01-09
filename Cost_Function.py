# ------------------------------------------------------------------------------
# The cost function is intentionally left in python script for general optimization problems
# as in case of modeling purpose for parameter tuning etc.
def Cost_Function(x, pDim):
    ##--------------------- Rosenbrock Banana Function---------------------------
    return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2

    ##-------------Shifted sphere function for multidimensional-------------------
    # cost = 0
    # for i in range(0,pDim):
    #      cost = cost + (x[i]-i)**2
    # return cost
