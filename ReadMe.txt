Use AISA optimization method with python and C-Extension:

--- Open your folder in terminal then run the following installation code to set up C extension
(C-extension is written to perform faster implementation of optimization code for large number of the parameters.) 

sudo CC=gcc python3 setup.py install

--- use run_AISA_cPython for the simple Rosenbrock function optimization
You can see the plot of cost function by removing the comments.

--- Change your problem definition (Cost function) and modify the Cost_Function function and 
initial parameters run_AISA_cPython file

--- Original paper ---

Bogar, Esref, and Selami Beyhan. "Adolescent Identity Search Algorithm
(AISA): A novel metaheuristic approach for solving optimization problems." Applied Soft 
Computing 95 (2020): 106503.

Remember that there is no guarantee to find the global minimum.

Good Luck!


