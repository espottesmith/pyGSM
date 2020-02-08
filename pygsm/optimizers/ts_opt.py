from pymatgen.core.structure import Molecule as PMGMol

from level_of_theories.xtb import XTB
from wrappers.molecule import Molecule
from potential_energy_surfaces.pes import PES
from .eigenvector_follow import eigenvector_follow
from utilities.manage_xyz import read_xyz, get_atoms, xyz_to_np
from coordinate_systems import Topology,PrimitiveInternalCoordinates,DelocalizedInternalCoordinates
from utilities import elements


def ts_opt(filepath, charge):

    pmg_mol = PMGMol.from_file(filepath)
    pmg_mol.set_charge_and_spin(charge=charge)

    # good
    geoms = read_xyz(filepath)
    lot = XTB.from_options(charge=charge, geom=geoms, nproc=32, ID=0)
    pes = PES.from_options(lot=lot)

    atom_symbols  = get_atoms(geoms)
    ELEMENT_TABLE = elements.ElementData()
    atoms = [ELEMENT_TABLE.from_symbol(atom) for atom in atom_symbols]
    xyz = xyz_to_np(geoms)
    top = Topology.build_topology(xyz,
                                  atoms,
                                  hybrid_indices=None,
                                  prim_idx_start_stop=None)

    p = PrimitiveInternalCoordinates.from_options(xyz=xyz, atoms=atoms, addtr=False,
                                                  addcart=False, topology=top)

    coord_obj = DelocalizedInternalCoordinates.from_options(xyz=xyz, atoms=atoms,
                                                            addtr=False, addcart=False,
                                                            primitives=p)

    print(type(coord_obj))

    mol = Molecule.from_options(geom=geoms, PES=pes, coord_obj=coord_obj,
                                Form_Hessian=True)

    optimizer = eigenvector_follow.from_options(
        print_level=2,
        Linesearch="NoLineSearch",
        update_hess_in_bg=False,
        conv_Ediff=100,
        conv_dE=1,
        conv_gmax=100,
        opt_climb=False,
        )

    geoms, energies = optimizer.optimize(mol,
                                         refE=mol.energy,
                                         opt_type="TS",
                                         opt_steps=1000)

    return geoms, energies

