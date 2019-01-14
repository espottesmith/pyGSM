import numpy as np
import options
import os
from se_gsm import *
from dlc import *
from penalty_pes import *
import pybel as pb
import sys


class SE_Cross(SE_GSM):

    @staticmethod
    def from_options(**kwargs):
        return SE_Cross(SE_Cross.default_options().set_values(kwargs))

    def go_gsm(self,max_iters=50,opt_steps=3):
        self.icoords[0].gradrms=0.
        self.icoords[0].energy = self.icoords[0].V0 = self.icoords[0].PES.get_energy(self.icoords[0].geom)
        print 'initial energy is {:1.4f}'.format(self.icoords[0].energy)
        sys.stdout.flush()
        self.interpolate(1)
        self.growth_iters(iters=max_iters,maxopt=opt_steps,nconstraints=self.nconstraints)

        print 'SE_Cross growth phase over'
        print 'Warning last node still not fully optimized'

        #if self.check_if_grown():
        #    self.icoords[self.nR] = DLC.copy_node(self.icoords[self.nR-1],self.nR,0)
        #    self.nR += 1

        if self.tstype==3:
            self.icoords[self.nR-1].OPTTHRESH=0.01
            oiters=50
            nconstraints=1
            ictan = DLC.tangent_SE(self.icoords[self.nR-1],self.driving_coords)
            #ictan = DLC.tangent_1(self.icoords[self.nR-1],self.icoords[self.nR-2])
        else:
            self.icoords[self.nR-1].OPTTHRES=self.CONV_TOL
            oiters=100
            self.icoords[self.nR-1].opt_type=0
            nconstraints=0
            ictan=None

        self.optimize(n=self.nR-1,nsteps=oiters,nconstraints=nconstraints,ictan=ictan)
            
        if self.tstype==3:
            self.write_xyz_files(iters=1,base="after_penalty",nconstraints=self.nconstraints)
            self.icoords[self.nR] = DLC.copy_node_X(self.icoords[self.nR-1],new_node_id=self.nR,rtype=5)
            self.icoords[self.nR].OPTTHRESH=self.CONV_TOL
            self.optimize(n=self.nR,nsteps=oiters*2,nconstraints=2)
            self.write_xyz_files(iters=1,base="grown_string",nconstraints=self.nconstraints)
        else:
            self.write_xyz_files(iters=1,base="grown_string",nconstraints=self.nconstraints)

    def add_node(self,n1,n2,n3=None):
        print "adding node: %i from node %i"%(n2,n1)
        return DLC.add_node_SE_X(self.icoords[n1],self.driving_coords,dqmag_max=self.DQMAG_MAX,dqmag_min=self.DQMAG_MIN)
    
    def converged(self,n):

        if self.icoords[n].opt_type==0:
            tmp1 = np.copy(self.icoords[n].PES.grad1)
            tmp2 = np.copy(self.icoords[n].PES.grad2)
            print 'norm1: {:1.4f} norm2: {:1.4f}'.format(np.linalg.norm(tmp1),np.linalg.norm(tmp2)),
            print 'ratio: {:1.4f}'.format(np.linalg.norm(tmp1)/np.linalg.norm(tmp2))
            tmp1 = tmp1/np.linalg.norm(tmp1)
            tmp2 = tmp2/np.linalg.norm(tmp2)
            print 'normalized gradient dot product:',float(np.dot(tmp1.T,tmp2))
            sys.stdout.flush()
            if self.icoords[n].gradrms<self.CONV_TOL and 1.-abs(float(np.dot(tmp1.T,tmp2))) <= 0.02 and abs(self.icoords[n].PES.dE) <= 1.25:
                return True
            else:
                return False
        elif self.icoords[n].opt_type==1: #constrained growth
            if self.icoords[n].gradrms<self.icoords[n].OPTTHRESH:
                return True
            else:
                return False
        elif self.icoords[n].opt_type==5:
            if self.icoords[n].gradrms<self.CONV_TOL and abs(self.icoords[n].PES.dE) <= 1.0:
                return True
            else:
                return False

    def check_if_grown(self):
        isDone = False
        epsilon = 1.5
        pes1dE = self.icoords[self.nR-1].PES.dE
        pes2dE = self.icoords[self.nR-2].PES.dE
        condition1 = (abs(self.icoords[self.nR-1].bdist) <=0.5* abs(self.icoords[1].bdist) and (abs(pes1dE) > abs(pes2dE)))
        condition2 = self.icoords[self.nR-1].bdist <=self.BDIST_MIN*(1.0-self.icoords[1].bdist)
        if condition1 or condition2:
            isDone = True
        return isDone

