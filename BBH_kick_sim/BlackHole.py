import numpy as np
import sympy
import surfinBH

class BlackHole:
    '''The Black hole object, used to keep track of their parameters'''
    def __init__(self,mass,spin,velocity,time,bh1=None,bh2=None,id=None):
        '''Creating a new Black hole object'''
        
        # Need to keep track of formation parameters and time evolved
        self.m = mass #Solar masses
        
        self.s = spin #Unitless spin
        self.s_mag = np.sqrt(np.sum(np.square(spin)))
        #self.s_i = spin #Does spin evolve over time?
        
        self.v = velocity #Km/s
        self.v_i = velocity
        
        self.t = time #s
        self.t_i = time
        
        self.parents = (bh1,bh2)
        
        self.id=id
    
    def n_parents(self):
        n = 0
        if self.parents[0]==None:
            return n
        n += 1 + self.parents[0].n_parents()
        n += 1 + self.parents[1].n_parents()
        
        return n
    
    def __str__(self):
        BH = f'[M:{round(self.m,3)},S_mag:{round(self.s_mag,3)},'
        BH+= f'v_mag:{round(np.sqrt(np.sum(np.square(self.v))),3)},'
        BH+= f'v_imag:{round(np.sqrt(np.sum(np.square(self.v_i))),3)},'
        BH+= f't:{round(self.t,3)},t_i:{round(self.t_i,3)}]'
        return BH
    
    def printTree(self,params=['m','s_mag']):
        from treelib import Node, Tree
        
        def _makeNode(tree,BH,final):
            parents = BH.parents
            
            parstring = '['
            if 'm' in params:
                parstring+=f'm={round(BH.m,3)},'
            if 's' in params:
                parstring+=f's={np.around(BH.s,3)},'
            if 's_mag' in params:
                parstring+=f's_mag={round(np.sqrt(np.sum(np.square(BH.v))),3)},'
            if 'v' in params:
                parstring+=f'v={round(BH.m,3)},'
            if 'v_mag' in params:
                parstring+=f'v_mag={round(np.sqrt(np.sum(np.square(BH.v))),3)},'
            if 'v_i' in params:
                parstring+=f'v_i={round(BH.v_i,3)},'
            if 'v_imag' in params:
                parstring+=f'v_imag={round(np.sqrt(np.sum(np.square(BH.v_i))),3)},'
            if 't' in params:
                parstring+=f't={BH.t},'
            if 't_i' in params:
                parstring+=f't_i={BH.t_i},'
            if 'id' in params:
                parstring+=f'id={BH.id},'
            parstring = parstring[:-1]+']'
            
            #Make Current node, then add nodes for each
            if final!=None:
                tree.create_node(parstring,BH.id,parent=final,data=BH)
            else:
                tree.create_node(parstring,BH.id,parent=None,data=BH)
            if parents[0]!=None:
                _makeNode(tree,parents[0],BH.id)
                _makeNode(tree,parents[1],BH.id)
        
        tree = Tree()
        _makeNode(tree,self,None)
        tree.show(key=lambda node: -node.data.m)

        
