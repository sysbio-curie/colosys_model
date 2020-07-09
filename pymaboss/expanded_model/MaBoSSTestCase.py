#!/usr/bin/env python
# coding: utf-8

# In[2]:


import unittest
import maboss
import math


# In[3]:


class MaBoSSTestCase(unittest.TestCase):
    
    def __init__(self, sim, verbose = True):
        unittest.TestCase.__init__(self)
        self.Simulation = sim.copy()
        self.Old_sim = (self.Simulation).copy()
        self.Old_result = None 
        self.New_sim = (self.Simulation).copy()     
        self.New_result = None
        self.VERBOSE = verbose
    
    #call at the end of the tests
    def resetSimulations(self):
        self.Old_sim = (self.Simulation).copy()
        self.Old_result = None
        self.New_sim = (self.Simulation).copy()
        self.New_result = None
    
    
    def checkNodes(self, state):
        for node in state:
            if node not in self.Simulation.network.keys():
                print(node, 'is not present in the network.')
                return False
        
        return True
    
    
    #set output for save time when you run a simulation   
    def setOutput(self, output):
        self.Old_sim.network.set_output(output)
        self.New_sim.network.set_output(output)
              
        
    def setInitialConditions(self, initial_conditions):
        
        if not initial_conditions and self.VERBOSE: 
            print('Care, you did not changed any initial condition!')
            return
        
        if type(initial_conditions) is dict:
            for node in initial_conditions:
                self.Old_sim.network.set_istate(node, initial_conditions[node])
                self.New_sim.network.set_istate(node, initial_conditions[node])
        #input: list: list[0]=list of node, list[1] dictionary of prob        
        elif type(initial_conditions) is list:
            self.Old_sim.network.set_istate(initial_conditions[0], initial_conditions[1])
            self.New_sim.network.set_istate(initial_conditions[0], initial_conditions[1])
        
        else: print('Invalid initial condition!')
            
        return
    
    
    #mutate New_simulation if a list of mutations is given
    def mutateSimulation(self, mutations={}):
        if mutations: 
            for n in mutations:
                (self.New_sim).mutate(n,mutations[n])
                
        elif self.VERBOSE: print('Care, you did not insert mutations!')
            
        
    def runBothSimulations(self):
        self.Old_result = self.Old_sim.run()
        self.New_result = self.New_sim.run()
        
    #return the states of a sim with their probability
    def getProbtrajStates(self, sim_name):
        if sim_name == 'New': result = (self.New_result).get_last_states_probtraj()
        elif sim_name == 'Old': result = (self.Old_result).get_last_states_probtraj()
        
        states = result.columns    
        probability_states = {}  # dict----> key:probability  value:state in list form
    
        for state in states:
            prob = result[state].values[0]
            probability_states[prob] = state.split(' -- ')
            
        return probability_states
    
    #return the stable states of a sim with their probability
    def getStableStates(self, sim_name):
        if sim_name == 'New': s_states_table = self.New_result.get_fptable()
        elif sim_name == 'Old': s_states_table = self.Old_result.get_fptable()
            
        s_states = s_states_table['State']
        s_states = [state.split(' -- ') for state in s_states]
        s_states_prob = s_states_table['Proba']
        
        probability_states = dict(zip(s_states_prob,s_states))
        
        return probability_states
            
    #return a dictionari with the states satisfing the condition with probability as key     
    def checkForState(self, kind, condition={}, all_states={}):
        
        #check if I'm looking for <nil> state only for last prob traj :CHECK CASE IN WICH ALL STATES ARE 0
        if ( kind == 'last' and all(i == 0 for i in condition.values()) ): 
            active_nodes = ['<nil>']
            inactive_nodes = []
        
        else:
            active_nodes = [x for x in condition if condition[x]==1]
            inactive_nodes = [x for x in condition if condition[x]==0]
        
        resulting_states = all_states.copy()    #states that s
        
        for node in active_nodes:
            resulting_states = {x: resulting_states[x] for x in resulting_states if node in resulting_states[x]}
            
        for node in inactive_nodes:
            resulting_states = {x: resulting_states[x] for x in resulting_states if node not in resulting_states[x]}
        
        if resulting_states: return resulting_states  #if there is not match the probability is 0
        else: return None #{0: condition}
        
        
            
    def truncate(self, number, digits) -> float:
        stepper = 10.0 ** digits
        return math.trunc(stepper * number) / stepper
    
    #print the dictionary: prob_ state in a nice form
    def printStates(self, states):
        if states == None: print('No one')
        else:
            for prob in states:
                print('\nProbability = ', prob, '\nState: ', str(states[prob]))
        
        return
    
    
    def compare(self, Old_p, New_p, direction, message, digits):
        Old_p = self.truncate(Old_p, digits)
        New_p = self.truncate(New_p, digits)
              
        #INCREASE?
        if direction == self.INCREASE: 
            try:
                self.assertTrue(New_p > Old_p)
                if self.VERBOSE: print("True! ", message)
                else: print('... OK')
            
            except Exception as e:
                if self.VERBOSE: print("False! ", message)
                else: print(e)
                
        #DECREASE?   
        elif direction == self.DECREASE:
            try:
                self.assertTrue(New_p < Old_p)
                if self.VERBOSE: print("True! ", message)
                else: print('... OK')
            
            except Exception as e:
                if self.VERBOSE: print("False! ", message)
                else: print(e)
                
        #STAY THE SAME?
        elif direction == self.CHANGE: 
            try:
                self.assertAlmostEqual(Old_p, New_p)
                if self.VERBOSE: print("True! ", message)
                else: print('... OK')
            
            except Exception as e:
                if self.VERBOSE: print("False! ", message)
                else: print(e)
            
           
        return   
    
    
    ##################################################################################################################
    ############################################ ASSERTION FUNCTIONS #################################################
    ##################################################################################################################
    #take last nodes probability for every node
    def getFullLastNodesProbabilities(self, mutations, I_C):
        
        #print('The simulation will be run with all the nodes external')
        
        output = list(self.New_sim.network.keys())
        
        self.setOutput(output)
        self.setInitialConditions(I_C)
        self.mutateSimulation(mutations)
        
        self.New_result = self.New_sim.run()
        
        prob_table = self.New_result.get_last_nodes_probtraj()
        
        
        probability = {}
    
        for node in prob_table.columns:
            p = prob_table[node].values[0]
            p = self.truncate(p,10)
            probability[node]=[1-p,p]
            
        
        for node in self.New_sim.network.keys():
            if node not in probability.keys():
                probability[node]=[1,0]
                
                
        self.resetSimulations()
        
        return probability
    
    def assertStateProbability(self, mutations, I_C, state, direction, digits = 4):
        
        if not self.checkNodes(state): return
        
        output = list(state.keys())
        
        self.setOutput(output)             #set the output for Old/New_sim 
        self.setInitialConditions(I_C)     #set initial conditions
        self.mutateSimulation(mutations)   #mutate New_sim
        
        self.runBothSimulations()
        
        all_Old_states = self.getProbtrajStates('Old')  #dict: key:stateProb  value:state (list form)
        all_New_states = self.getProbtrajStates('New')
        
        Old_states = self.checkForState('last', state, all_Old_states)  #there should be prob_state:state_list_form 
        New_states = self.checkForState('last', state, all_New_states)  #if prob is not 0  
        
        
        if Old_states == None: Old_state_prob = 0
        else: Old_state_prob = list(Old_states.keys())[0]                 #if possible invert keys and value
            
        if New_states == None: New_state_prob = 0
        else: New_state_prob = list(New_states.keys())[0]
        
        #check it is all fine
        if (Old_states != None and New_states != None and len(Old_states)!=1 and len(Old_states)!=1):
            print('ERROR, TO MUTCH STATES')
            return
         
        message = '\nThe new probability of reaching the state is: {}'         ' \nThe old one is: {}'.format(New_state_prob, Old_state_prob)
        
        self.compare(Old_state_prob, New_state_prob, direction, message, digits) #ad an output!!!!!
        
        self.resetSimulations()
    
        return
    
    
    def assertStableStateProbability(self, mutations, I_C, state, direction, digits = 4):
         
        if not self.checkNodes(state): return
        
        self.setInitialConditions(I_C)    #set initial conditions
        self.mutateSimulation(mutations)   #mutate New_sim
        
        self.runBothSimulations()
        
        all_Old_stable_states = self.getStableStates('Old')
        all_New_stable_states = self.getStableStates('New')
        
        Old_states = self.checkForState('stable', state, all_Old_stable_states)
        New_states = self.checkForState('stable', state, all_New_stable_states)     
        
        Old_state_prob = 0
        if Old_states != None: 
            for prob in Old_states.keys():
                Old_state_prob += prob
        
        New_state_prob = 0
        if New_states != None: 
            for prob in New_states.keys():
                New_state_prob += prob
            
        message = '\nThe new probability of reaching the stable state is: {}'         ' \nThe old one is: {}'.format(New_state_prob, Old_state_prob)
        
        self.compare(Old_state_prob, New_state_prob, direction, message, digits)
            
        self.resetSimulations()
        
        return
    
    
    #NAME?????
    #maBoss only find ss with prob >0?
    #assert if given some condition (nodes dictionaries) other nodes has always the same values
    def assertStableStates(self, mutations, condition, nodes_expected):
        
        if not self.checkNodes(condition): return
        if not self.checkNodes(nodes_expected): return
        
        #for stable states don't need outputs
        self.mutateSimulation(mutations)   #mutate New_sim
        
        self.New_result = self.New_sim.run()
        
        all_stable_states = self.getStableStates('New')
        #states satisfying the input condition
        states_satisfying_condition = self.checkForState('stable', condition, all_stable_states)
        #selected states satisfying both condition and nodes_expected
        if states_satisfying_condition != None:
            if not nodes_expected: selected_states = None 
            else: selected_states = self.checkForState('stable', nodes_expected, states_satisfying_condition)
        
        else:
            print("Not even one stable state satisfy: ", condition)
            return
        
        try:
            self.assertEqual(states_satisfying_condition, selected_states)
            if self.VERBOSE: print("True! \nAll the states that satisfy: ", condition, " have: ", nodes_expected)
            else: print("... OK")
                
        except Exception as e: 
            if self.VERBOSE: 
                print('False! \nThe states with ', condition, ' are : ') 
                self.printStates(states_satisfying_condition)
                print('\nOf these, those with ', nodes_expected, ' are: ')                        
                self.printStates(selected_states)                  
                                         
            else: print(e)
            
        
        self.resetSimulations()
        
        return
    
    
    #IT IS CLEAR??
    #compare the probability of reaching a state with prob (1 ask INCREASE, -1 DECREASE, 0 IS EQUAL)
    def compareStableStateProbability(self, mutations, state, direction, reference_prob = 0, digits = 4):
         
        if not self.checkNodes(state): return    
        
        self.mutateSimulation(mutations)          #mutate New_sim
        
        self.New_result = self.New_sim.run()
        
        all_New_stable_states = self.getStableStates('New')
        
        New_states = self.checkForState('stable', state, all_New_stable_states)     
        
        New_state_prob = 0
        if New_states != None: 
            for prob in New_states.keys():
                New_state_prob += prob
            
        message = '\nThe probability of reaching the state is: {}'.format(New_state_prob)
        print('The reference probability is: ', reference_prob)
        self.compare(reference_prob, New_state_prob, direction, message, digits)
            
        self.resetSimulations()
        
        return
        
        
        
    INCREASE = 'increase'
    DECREASE = 'decrease'
    CHANGE = 'stable'                      #actually it is: does not change  
    


# In[ ]:




