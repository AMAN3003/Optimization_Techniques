#!/usr/bin/env python
#
# -*- Mode: python -*-

""" Simplex - a regression method for arbitrary nonlinear function minimization

Simplex minimizes an arbitrary nonlinear function of N variables by the
Nedler-Mead Simplex method as described in:

It makes no assumptions about the smoothness of the function being minimized.
It converges to a local minimum which may or may not be the global minimum
depending on the initial guess used as a starting point.
"""

import math
import copy

class Simplex:
    def __init__(self, func_tester, guess_list, increment_list, ref_const = -1, eqx_const = 2, contract_const = 0.5):
        """Simplex is initialized.
        INPUTS
        ------
        func_tester       minimize function
        guess_list[]       initial guess list
        increment_list[]  increments list with ,  perturbation size
        ref_const            reflection constant used in simplex
        eqx_const            expansion constant used in simplex
        contract_const       contraction constant used in simplex
        """
        self.func_tester = func_tester
        self.guess_list = guess_list
        self.increment_list = increment_list
        self.ref_const = ref_const
        self.eqx_const = eqx_const
        self.contract_const = contract_const
        self.vars_num = len(self.guess_list)


    def minimize(self, simplex_epsilon = 0.0001, Iter_maxnum = 250, Output_solved = 1, **kwargs):
        """Walks through the simplex down to a local minima.at which it converges
        INPUTS
        ------
        simplex_epsilon       for the convergence its required
        Iter_maxnum           for the maximum number of iterations allowed
        Output_solved     	  check if has optimal solution that is output is there  

        OUTPUTS
        -------
        
        an array with the final value,lowest value of the simplex error funtion and the no of iteration to get there

        """
        self.simplex_list = []

        self.lowest_val = -1
        self.highest_val = -1
        self.Secondhigh_val = -1

        self.List_Errors = []
        self.Errorval_Present = 0
        # Initialization of the vertixes 
        for Simplex_Vertex in range(0, self.vars_num + 3): # 2 are extra for centroid and reflexion point
            self.simplex_list.append(copy.copy(self.guess_list))
        # increment for the list is used 
        for Simplex_Vertex in range(0, self.vars_num + 1):
            for x in range(0, self.vars_num):
                if x == (Simplex_Vertex - 1):
                    self.simplex_list[Simplex_Vertex][x] = self.guess_list[x] + self.increment_list[x]
            self.List_Errors.append(0)
        self.Error_vertex_Cal(**kwargs)

        iter = 0
        for iter in range(0, Iter_maxnum):
            # Identify highest, second highest, and lowest vertices
            self.highest_val = 0
            self.lowest_val = 0
            for Simplex_Vertex in range(0, self.vars_num + 1):
                if self.List_Errors[Simplex_Vertex] > self.List_Errors[self.highest_val]:
                    self.highest_val = Simplex_Vertex
                if self.List_Errors[Simplex_Vertex] < self.List_Errors[self.lowest_val]:
                    self.lowest_val = Simplex_Vertex
            self.Secondhigh_val = 0
            for Simplex_Vertex in range(0, self.vars_num + 1):
                if Simplex_Vertex == self.highest_val:
                    continue
                if self.List_Errors[Simplex_Vertex] > self.List_Errors[self.Secondhigh_val]:
                    self.Secondhigh_val = Simplex_Vertex
            #convergence test implemented 
            Simplexvar = 0.0
            Simplexvar1 = 0.0
            for Simplex_Vertex in range(0, self.vars_num + 1):
                Simplexvar = Simplexvar + self.List_Errors[Simplex_Vertex]
            Funcvar = Simplexvar / (self.vars_num + 1)
            for Simplex_Vertex in range(0, self.vars_num + 1):
                Simplexvar1 = Simplexvar1 + (self.List_Errors[Simplex_Vertex] - Funcvar)**2
            testval = math.sqrt(Simplexvar1 / self.vars_num)
            
            # it print the progress of the solution
            if Output_solved:
		    print '#%d: Bestsol = %f   Worstsol = %f' % (iter,self.List_Errors[self.lowest_val],self.List_Errors[self.highest_val]),
		    for Simplex_Vertex in range(0, self.vars_num + 1):
			    print "[",
			    for x in range(0, self.vars_num):
				print "%.2f" % self.simplex_list[Simplex_Vertex][x],
			    print "]",
		    print

                
            if testval <= simplex_epsilon:   # checked converged  so loop break
                break;
            else:                   # Didn't converge. keep going on
                # Calculate centroid for the simplex,ignoring the highest simplex vertex 
                for x in range(0, self.vars_num):
                    Simplexvar = 0.0
                    for Simplex_Vertex in range(0, self.vars_num + 1):
                        if Simplex_Vertex == self.highest_val:
                            continue
                        Simplexvar = Simplexvar + self.simplex_list[Simplex_Vertex][x]
                    self.simplex_list[self.vars_num + 1][x] = Simplexvar / self.vars_num

                self.simplex_reflection_funct()

                self.Errorval_Present = self.func_tester(self.guess_list, **kwargs)

                if self.Errorval_Present < self.List_Errors[self.lowest_val]:
                    tmpval = self.Errorval_Present
                    self.simplex_Expand_funct()
                    self.Errorval_Present = self.func_tester(self.guess_list, **kwargs)
                    if self.Errorval_Present < tmpval:
                        self.Right_Expansion_Point()
                    else:
                        self.Errorval_Present = tmpval
                        self.Right_Reflexion_Point()

                elif self.Errorval_Present <= self.List_Errors[self.Secondhigh_val]:
                    self.Right_Reflexion_Point()

                elif self.Errorval_Present <= self.List_Errors[self.highest_val]:
                    self.Right_Reflexion_Point()

                    self.simplex_Contact_funct()
                    self.Errorval_Present = self.func_tester(self.guess_list, **kwargs)
                    if self.Errorval_Present < self.List_Errors[self.highest_val]:
                        self.Right_Contaction_Point()
                    else:
                        self.simplex_Contraction_left(**kwargs)

                elif self.Errorval_Present >= self.List_Errors[self.highest_val]:
                    self.simplex_Contact_funct()
                    self.Errorval_Present = self.func_tester(self.guess_list, **kwargs)
                    if self.Errorval_Present < self.List_Errors[self.highest_val]:
                        self.Right_Contaction_Point()
                    else:
                        self.simplex_Contraction_left(**kwargs)

        # case when converged or maximum no of iteration reached
        # this returns  lowest simplex vertex and the presenterrors.
        for x in range(0, self.vars_num):
            self.guess_list[x] = self.simplex_list[self.lowest_val][x]
        self.Errorval_Present = self.List_Errors[self.lowest_val]
        return self.guess_list, self.Errorval_Present, iter

    def simplex_Contact_funct(self):
        for x in range(0, self.vars_num):
            self.guess_list[x] = self.contract_const * self.simplex_list[self.highest_val][x] + (1 - self.contract_const) * self.simplex_list[self.vars_num + 1][x]
        return

    def simplex_Expand_funct(self):
        for x in range(0, self.vars_num):
            self.guess_list[x] = self.eqx_const * self.guess_list[x]                 + (1 - self.eqx_const) * self.simplex_list[self.vars_num + 1][x]
        return

    def simplex_reflection_funct(self):
        for x in range(0, self.vars_num):
            self.guess_list[x] = self.ref_const * self.simplex_list[self.highest_val][x] + (1 - self.ref_const) * self.simplex_list[self.vars_num + 1][x]
            self.simplex_list[self.vars_num + 2][x] = self.guess_list[x] # REMEMBER THE REFLECTED POINT
        return

    def simplex_Contraction_left(self, **kwargs):
        for Simplex_Vertex in range(0, self.vars_num + 1):
            if Simplex_Vertex == self.lowest_val:
                continue
            for x in range(0, self.vars_num):
                self.simplex_list[Simplex_Vertex][x] = 0.5 * (self.simplex_list[Simplex_Vertex][x] + self.simplex_list[self.lowest_val][x])
        self.Error_vertex_Cal(**kwargs)
        return

    def Right_Contaction_Point(self):
        self.List_Errors[self.highest_val] = self.Errorval_Present
        for x in range(0, self.vars_num):
            self.simplex_list[self.highest_val][x] = self.guess_list[x]
        return

    def Right_Expansion_Point(self):
        self.List_Errors[self.highest_val] = self.Errorval_Present
        for x in range(0, self.vars_num):
            self.simplex_list[self.highest_val][x] = self.guess_list[x]
        return

    def Right_Reflexion_Point(self):
        self.List_Errors[self.highest_val] = self.Errorval_Present
        for x in range(0, self.vars_num):
            self.simplex_list[self.highest_val][x] = self.simplex_list[self.vars_num + 2][x]
        return

    def Error_vertex_Cal(self,**kwargs):
        for Simplex_Vertex in range(0, self.vars_num + 1):
            if Simplex_Vertex == self.lowest_val:
                continue
            for x in range(0, self.vars_num):
                self.guess_list[x] = self.simplex_list[Simplex_Vertex][x]
            self.Errorval_Present = self.func_tester(self.guess_list, **kwargs)
            self.List_Errors[Simplex_Vertex] = self.Errorval_Present
        return

def Simplex_objective_function(args, static_par):
	return abs(- (args[0]-static_par[0])**2 - (args[1]-static_par[1])**2)
	
def main():
	Simplexobj = Simplex(Simplex_objective_function, [1, 1, 1], [.01, .01, .01])
	values, err, iter = Simplexobj.minimize(simplex_epsilon = 0.0000001, Iter_maxnum = 250, Output_solved = 0, static_par=[2,4])
	del(Simplexobj)
	print 'arguments = ', values
	print 'errors = ', err
	print 'no of iterations = ', iter
	
if __name__ == '__main__':
    main()

