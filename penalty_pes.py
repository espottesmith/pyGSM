import options
from pes import * 
import sys

class Penalty_PES(PES):
    """ penalty potential energy surface calculators """

    def __init__(self,
            PES1,
            PES2,
            sigma=3.5,
            alpha=0.02*KCAL_MOL_PER_AU):
        self.PES1 = PES1
        self.PES2 = PES2
        self.lot = PES1.lot
        self.alpha = alpha
        self.sigma = sigma
        print 'PES1 multiplicity: {} PES2 multiplicity: {}'.format(self.PES1.multiplicity,self.PES2.multiplicity)


    @staticmethod
    def from_options(**kwargs):
        """ Returns an instance of this class with default options updated from values in kwargs"""
        return Penalty_PES(Penalty_PES.default_options().set_values(kwargs))

    def get_energy(self,geom):
        E1 = self.PES1.get_energy(geom)
        E2 = self.PES2.get_energy(geom)
        #avgE = 0.5*(self.PES1.get_energy(geom) + self.PES2.get_energy(geom))
        avgE = 0.5*(E1+E2)
        #self.dE = self.PES2.get_energy(geom) - self.PES1.get_energy(geom)
        self.dE = E2-E1
        print "E1: %1.4f E2: %1.4f"%(E1,E2),
        print "delta E = %1.4f" %self.dE,
        #TODO what to do if PES2 is or goes lower than PES1?
        G = (self.dE*self.dE)/(abs(self.dE) + self.alpha)
        #if self.dE < 0:
        #    G*=-1
        print "G = %1.4f" % G
        print "alpha: %1.4f sigma: %1.4f"%(self.alpha,self.sigma),
        print "F: %1.4f"%(avgE+self.sigma*G)
        sys.stdout.flush()
        return avgE+self.sigma*G

    def get_gradient(self,geom):
        self.grad1 = self.PES1.get_gradient(geom)
        self.grad2 = self.PES2.get_gradient(geom)
        avg_grad = 0.5*(self.grad1 + self.grad2)
        dgrad = self.grad2 - self.grad1
        if self.dE < 0:
            dgrad *= -1
        factor = self.sigma*((self.dE*self.dE) + 2.*self.alpha*abs(self.dE))/((abs(self.dE) + self.alpha)**2)
        #print "factor is %1.4f" % factor
        return avg_grad + factor*dgrad


if __name__ == '__main__':
    if 1:
        from pytc import *
        import dlc as ic
        import pybel as pb
        from units import  *
        #filepath1="tests/ethylene.xyz"
        filepath2="tests/twisted_ethene.xyz"
        nocc=7
        nactive=2
        lot=PyTC.from_options(states=[(1,0),(1,1)],nocc=nocc,nactive=nactive,basis='6-31gs')
        lot.cas_from_file(filepath2)
        #lot.casci_from_file_from_template(filepath1,filepath2,nocc,nocc)
        pes1 = PES.from_options(lot=lot,ad_idx=0,multiplicity=1)
        pes2 = PES.from_options(lot=lot,ad_idx=1,multiplicity=1)

        p = Penalty_PES(pes1.options.copy().set_values({
            'PES1':pes1,
            'PES2':pes2,
            }))
        #p.get_energy()
        #p.getGrad()
        mol=pb.readfile("xyz",filepath2).next()
        mol.OBMol.AddBond(6,1,1)
        print "ic1"
        ic1=ic.ICoord.from_options(mol=mol,PES=p,resetopt=False)

        for i in range(1):
            ic1.optimize(100)
            if ic1.gradrms<ic1.OPTTHRESH:
                break
            if  ic1.lot.dE>0.001*KCAL_MOL_PER_AU:
                ic1.lot.sigma *=2.
                print "increasing sigma %1.2f" % ic1.lot.sigma
        print "Finished"
    if 0:
        from molpro import *
        import icoord as ic
        import pybel as pb
        from units import  *
        filepath="tests/twisted_ethene.xyz"
        geom=manage_xyz.read_xyz(filepath,scale=1)   
        nocc=7
        nactive=2
        lot=Molpro.from_options(states=[(1,0),(1,1)],nocc=nocc,nactive=nactive,basis='6-31g*',nproc=4)
        pes1 = PES.from_options(lot=lot,ad_idx=0,multiplicity=1)
        pes2 = PES.from_options(lot=lot,ad_idx=1,multiplicity=1)
        #p = Penalty_PES(pes1,pes2)
        p = Penalty_PES(pes1.options.copy().set_values({
            'PES1':pes1,
            'PES2':pes2,
            }))

        #e =p.get_energy(geom)
        #g =p.get_gradient(geom)
        #print e 
        #print g
        mol=pb.readfile("xyz",filepath).next()
        mol.OBMol.AddBond(6,1,1)
        print "####### ic1 ##########"
        ic1=ic.ICoord.from_options(mol=mol,PES=p,resetopt=False)

