# standard library imports
import sys
import os
from os import path

# third party
import numpy as np

# local application imports
sys.path.append(path.dirname( path.dirname( path.abspath(__file__))))
from .base_lot import Lot
from utilities import *


class XTB(Lot):
    def __init__(self,options):
        super(XTB,self).__init__(options)

        scratch = os.environ("SCRATCH")

        for state in self.states:
            tempfolder = scratch + '/string_{:03d}/{}.{}/'.format(self.ID,self.node_id,state[0])
            print(" making temp folder {}".format(tempfolder))
            os.system('mkdir -p {}'.format(tempfolder))

        copy_input_file = os.getcwd() + "/xtbinput.txt"
        print(copy_input_file)
        self.write_preamble(self.geom,self.states[0][0],copy_input_file)

    def write_preamble(self,geom,multiplicity,tempfilename):

        with open(tempfilename, "w") as tempfile:
            if not self.lot_inp_file:
                tempfile.write('$gfn\n')
                tempfile.write('\tmethod=2\n')
                tempfile.write('\tscc=true\n')
                tempfile.write('$scc\n')
                tempfile.write('\tmaxiterations=250\n')
                tempfile.write('\ttemp=300\n')
                tempfile.write('$write\n')
                tempfile.write('\tcharges=true\n')
                tempfile.write('\tmulliken=true\n')
                tempfile.write('\tdistances=true\n')
                tempfile.write('\tangles=true\n')
                tempfile.write('\ttorsions=true\n')
                tempfile.write('\tgbsa=true\n')
                tempfile.write('$end\n')
            else:
                with open(self.lot_inp_file) as lot_inp:
                    lot_inp_lines = lot_inp.readlines()
                for line in lot_inp_lines:
                    tempfile.write(line)

        with open(os.path.join(os.path.dirname(tempfilename), "coords"), 'w') as coordsfile:
            coordsfile.write("{}\n".format(len(geom)))
            coordsfile.write("{} {}\n").format(self.charge, multiplicity)

            if os.path.isfile("link.txt"):
                with open("link.txt") as link:
                    link_lines = link.readlines()
                tmp_geom = [list(i) for i in geom]
                for i,coord in enumerate(tmp_geom):
                    coord.append(link_lines[i].rstrip('\n'))
                    for i in coord:
                        coordsfile.write(str(i)+' ')
                    coordsfile.write('\n')
            else:
                for coord in geom:
                    for i in coord:
                        coordsfile.write(str(i)+' ')
                    coordsfile.write('\n')

    def run(self,geom,multiplicity):

        scratch = os.environ['SCRATCH']
        tempdir = scratch + '/string_{:03d}/{}.{}/'.format(self.ID,self.node_id,multiplicity)
        tempfilename = os.path.join(tempdir, 'tempxtbinp')

        self.write_preamble(geom,multiplicity,tempfilename)

        if self.calc_grad:
            cmd = "xtb {} --grad --gbsa water extreme --chrg {} --acc 0.01 --input {} --namespace {} --parallel 32".format(
                os.path.join(tempdir, "coords"),
                self.charge,
                tempfilename,
                os.path.join(tempdir, "xtbout")
            )
        else:
            cmd = "xtb {} --sp --gbsa water extreme --chrg {} --acc 0.01 --input {} --namespace {} --parallel 32".format(
                os.path.join(tempdir, "coords"),
                self.charge,
                tempfilename,
                os.path.join(tempdir, "xtbout")
            )

        #print(cmd)
        os.system(cmd)

        # PARSE OUTPUT #
        if self.calc_grad:
            epath = os.path.join(tempdir, "xtbout.energy")
            gradpath = os.path.join(tempdir, "xtbout.gradient")
            with open(epath) as efile:
                elines = efile.readlines()
            with open(gradpath) as gradfile:
                glines = gradfile.readlines()

            for line in elines:
                if "$" not in line:
                    self.E.append((multiplicity, float(line.split()[-1])))
                    break

            tmplist = glines[-(len(geom) + 1): -1]
            glist = list()
            for line in tmplist:
                split = [e.replace("D", "E") for e in line.split()]
                glist.append([float(i) for i in split])

            self.grada.append((multiplicity, glist))

        else:
            raise NotImplementedError

        return

    def get_energy(self,coords,multiplicity,state):
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = manage_xyz.np_to_xyz(self.geom,self.currentCoords)
            self.runall(geom)
            self.hasRanForCurrentCoords=True

        tmp = self.search_tuple(self.E,multiplicity)
        return np.asarray(tmp[state][1])*units.KCAL_MOL_PER_AU

    def get_gradient(self,coords,multiplicity,state):
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = manage_xyz.np_to_xyz(self.geom,self.currentCoords)
            self.runall(geom)
        tmp = self.search_tuple(self.grada,multiplicity)
        return np.asarray(tmp[state][1])*units.ANGSTROM_TO_AU

    @classmethod
    def copy(cls,lot,options,copy_wavefunction=True):
        base = os.environ['SCRATCH']
        node_id = options.get('node_id',1)

        if node_id != lot.node_id:  #and copy_wavefunction: # other theories are more sensitive than qchem -- commenting out
            for state in lot.states:
                multiplicity = state[0]
                efilepath_old=base+ '/string_{:03d}/{}.{}'.format(lot.ID,lot.node_id,multiplicity)
                efilepath_new =base+ '/string_{:03d}/{}.{}'.format(lot.ID,node_id,multiplicity)
                cmd = 'cp -r ' + efilepath_old +' ' + efilepath_new
                print(" copying SCRATCH files\n {}".format(cmd))
                os.system(cmd)
        return cls(lot.options.copy().set_values(options))

