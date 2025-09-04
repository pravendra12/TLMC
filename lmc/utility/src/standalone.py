import numpy as np
import spglib
from typing import Sequence, TypeVar
from ase import Atoms
from ase.data import chemical_symbols

from math import floor, log10

Vector = list[float]
T = TypeVar("T")


def get_scaled_positions(
    positions: np.ndarray,
    cell: np.ndarray,
    wrap: bool = True,
    pbc: list[bool] = [True, True, True],
) -> np.ndarray:
    """
    Returns the positions in reduced (or scaled) coordinates.

    Parameters
    ----------
    positions
        Atomic positions in Cartesian coordinates.
    cell
        Cell metric.
    wrap
        If ``True`` positions outside the unit cell will be wrapped into
        the cell in the directions with periodic boundary conditions
        such that the scaled coordinates are between zero and one.
    pbc
        Periodic boundary conditions.
    """

    fractional = np.linalg.solve(cell.T, positions.T).T

    if wrap:
        for i, periodic in enumerate(pbc):
            if periodic:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test.
                fractional[:, i] %= 1.0
                fractional[:, i] %= 1.0

    return fractional


def get_primitive_structure(
    structure: Atoms,
    no_idealize: bool = True,
    to_primitive: bool = True,
    symprec: float = 1e-5,
) -> Atoms:
    """
    Returns the primitive structure using spglib.

    Parameters
    ----------
    structure
        Input atomic structure.
    no_idealize
        If ``True`` lengths and angles are not idealized.
    to_primitive
        If ``True`` convert to primitive structure.
    symprec
        Tolerance imposed when analyzing the symmetry using spglib.
    """
    structure_copy = structure.copy()
    structure_as_tuple = ase_atoms_to_spglib_cell(structure_copy)

    ret = spglib.standardize_cell(
        structure_as_tuple,
        to_primitive=to_primitive,
        no_idealize=no_idealize,
        symprec=symprec,
    )
    if ret is None:
        raise ValueError(
            "spglib failed to find the primitive cell, maybe caused by large symprec."
        )
    lattice, scaled_positions, numbers = ret
    scaled_positions = [np.round(pos, 12) for pos in scaled_positions]
    structure_prim = Atoms(
        scaled_positions=scaled_positions,
        numbers=numbers,
        cell=lattice,
        pbc=structure.pbc,
    )
    structure_prim.wrap()

    return structure_prim


'''
def get_fractional_positions_from_neighbor_list(
        structure: Structure, neighbor_list: list) -> list[Vector]:
    """
    Returns the fractional positions of the lattice sites in structure from
    a neighbor list.

    Parameters
    ----------
    structure
        Input atomic structure.
    neighbor_list
        List of lattice neighbors of the input structure.
    """
    neighbor_positions = []
    fractional_positions = []
    lattice_site = LatticeSite(0, [0, 0, 0])

    for i in range(len(neighbor_list)):
        lattice_site.index = i
        position = structure.get_position(lattice_site)
        neighbor_positions.append(position)
        for neighbor in neighbor_list[i]:
            position = structure.get_position(neighbor)
            neighbor_positions.append(position)

    if len(neighbor_positions) > 0:
        fractional_positions = get_scaled_positions(
            np.array(neighbor_positions),
            structure.cell, wrap=False,
            pbc=structure.pbc)

    return fractional_positions


def get_position_from_lattice_site(structure: Atoms, lattice_site: LatticeSite):
    """
    Returns the corresponding position from the lattice site.

    Parameters
    ---------
    structure
        Input atomic structure.
    lattice_site
        Specific lattice site of the input structure.
    """
    return structure[lattice_site.index].position + \
        np.dot(lattice_site.unitcell_offset, structure.cell)
'''


def fractional_to_cartesian(
    structure: Atoms, frac_positions: list[Vector]
) -> np.ndarray:
    """
    Converts fractional positions into Cartesian positions.

    Parameters
    ----------
    structure
        Input atomic structure.
    frac_positions
        Fractional positions.
    """
    return np.dot(frac_positions, structure.cell)


def get_permutation(container: Sequence[T], permutation: list[int]) -> Sequence[T]:
    """
    Returns the permuted version of container.
    """
    if len(permutation) != len(container):
        raise RuntimeError(
            "Container and permutation not of same size."
            f"{len(container)} != {len(permutation)}"
        )
    if len(set(permutation)) != len(permutation):
        raise Exception
    return [container[s] for s in permutation]


def ase_atoms_to_spglib_cell(
    structure: Atoms,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns a tuple comprising three components, corresponding to cell
    metric, atomic positions, and atomic species.

    Parameters
    ----------
    structure
        Input atomic structure.
    """
    return (
        structure.cell,
        structure.get_scaled_positions(),
        structure.get_atomic_numbers(),
    )


def get_occupied_primitive_structure(
    structure: Atoms, allowed_species: list[list[str]], symprec: float
) -> tuple[Atoms, list[tuple[str, ...]]]:
    """Returns an occupied primitive structure with hydrogen on
    sublattice 1, Helium on sublattice 2 and so on

    Parameters
    ----------
    structure
        Input structure.
    allowed_species
        Chemical symbols that are allowed on each site.
    symprec
        Tolerance imposed when analyzing the symmetry using spglib.
    """
    if len(structure) != len(allowed_species):
        raise ValueError("structure and chemical symbols need to be the same size.")
    sorted_symbols = sorted({tuple(sorted(s)) for s in allowed_species})

    decorated_primitive = structure.copy()
    for i, sym in enumerate(allowed_species):
        sublattice = sorted_symbols.index(tuple(sorted(sym))) + 1
        decorated_primitive[i].symbol = chemical_symbols[sublattice]

    decorated_primitive = get_primitive_structure(decorated_primitive, symprec=symprec)
    decorated_primitive.wrap()
    primitive_chemical_symbols: list[tuple[str, ...]] = []
    for atom in decorated_primitive:
        sublattice = chemical_symbols.index(atom.symbol)
        primitive_chemical_symbols.append(sorted_symbols[sublattice - 1])

    for symbols in allowed_species:
        if tuple(sorted(symbols)) in primitive_chemical_symbols:
            index = primitive_chemical_symbols.index(tuple(sorted(symbols)))
            primitive_chemical_symbols[index] = symbols
    return decorated_primitive, primitive_chemical_symbols


def atomic_number_to_chemical_symbol(numbers: list[int]) -> list[str]:
    """Returns the chemical symbols equivalent to the input atomic
    numbers.

    Parameters
    ----------
    numbers
        Atomic numbers.
    """

    symbols = [chemical_symbols[number] for number in numbers]
    return symbols


def chemical_symbols_to_numbers(symbols: list[str]) -> list[int]:
    """Returns the atomic numbers equivalent to the input chemical
    symbols.

    Parameters
    ----------
    symbols
        Chemical symbols.
    """

    numbers = [chemical_symbols.index(symbols) for symbols in symbols]
    return numbers


def get_wyckoff_sites(
    structure: Atoms,
    map_occupations: list[list[str]] = None,
    symprec: float = 1e-5,
    include_representative_atom_index: bool = False,
) -> list[str]:
    """Returns the Wyckoff symbols of the input structure. The Wyckoff
    sites are of general interest for symmetry analysis but can be
    especially useful when setting up, e.g., a
    :class:`SiteOccupancyObserver
    <mchammer.observers.SiteOccupancyObserver>`.
    The Wyckoff labels can be conveniently attached as an array to the
    structure object as demonstrated in the examples section below.

    By default the occupation of the sites is part of the symmetry
    analysis. If a chemically disordered structure is provided this
    will usually reduce the symmetry substantially. If one is
    interested in the symmetry of the underlying structure one can
    control how occupations are handled. To this end, one can provide
    the :attr:`map_occupations` keyword argument. The latter must be a
    list, each entry of which is a list of species that should be
    treated as indistinguishable. As a shortcut, if *all* species
    should be treated as indistinguishable one can provide an empty
    list. Examples that illustrate the usage of the keyword are given
    below.

    Parameters
    ----------
    structure
        Input structure. Note that the occupation of the sites is
        included in the symmetry analysis.
    map_occupations
        Each sublist in this list specifies a group of chemical
        species that shall be treated as indistinguishable for the
        purpose of the symmetry analysis.
    symprec
        Tolerance imposed when analyzing the symmetry using spglib.
    include_representative_atom_index
        If True the index of the first atom in the structure that is
        representative of the Wyckoff site is included in the symbol.
        This is in particular useful in cases when there are multiple
        Wyckoff sites sites with the same Wyckoff letter.

    Examples
    --------
    Wyckoff sites of a hexagonal-close packed structure::

        >>> from ase.build import bulk
        >>> structure = bulk('Ti')
        >>> wyckoff_sites = get_wyckoff_sites(structure)
        >>> print(wyckoff_sites)
        ['2d', '2d']


    The Wyckoff labels can also be attached as an array to the
    structure, in which case the information is also included when
    storing the Atoms object::

        >>> from ase.io import write
        >>> structure.new_array('wyckoff_sites', wyckoff_sites, str)
        >>> write('structure.xyz', structure)

    The function can also be applied to supercells::

        >>> structure = bulk('GaAs', crystalstructure='zincblende', a=3.0).repeat(2)
        >>> wyckoff_sites = get_wyckoff_sites(structure)
        >>> print(wyckoff_sites)
        ['4a', '4c', '4a', '4c', '4a', '4c', '4a', '4c',
         '4a', '4c', '4a', '4c', '4a', '4c', '4a', '4c']

    Now assume that one is given a supercell of a (Ga,Al)As
    alloy. Applying the function directly yields much lower symmetry
    since the symmetry of the original structure is broken::

        >>> structure.set_chemical_symbols(
        ...        ['Ga', 'As', 'Al', 'As', 'Ga', 'As', 'Al', 'As',
        ...         'Ga', 'As', 'Ga', 'As', 'Al', 'As', 'Ga', 'As'])
        >>> print(get_wyckoff_sites(structure))
        ['8g', '8i', '4e', '8i', '8g', '8i', '2c', '8i',
         '2d', '8i', '8g', '8i', '4e', '8i', '8g', '8i']

    Since Ga and Al occupy the same sublattice, they should, however,
    be treated as indistinguishable for the purpose of the symmetry
    analysis, which can be achieved via the :attr:`map_occupations`
    keyword::

        >>> print(get_wyckoff_sites(structure,
        ...       map_occupations=[['Ga', 'Al'], ['As']]))
        ['4a', '4c', '4a', '4c', '4a', '4c', '4a', '4c',
         '4a', '4c', '4a', '4c', '4a', '4c', '4a', '4c']

    If occupations are to ignored entirely, one can simply provide an
    empty list. In the present case, this turns the zincblende lattice
    into a diamond lattice, on which case there is only one Wyckoff
    site::

        >>> print(get_wyckoff_sites(structure, map_occupations=[]))
        ['8a', '8a', '8a', '8a', '8a', '8a', '8a', '8a',
         '8a', '8a', '8a', '8a', '8a', '8a', '8a', '8a']
    """
    structure_copy = structure.copy()
    if map_occupations is not None:
        if len(map_occupations) > 0:
            new_symbols = []
            for symb in structure_copy.get_chemical_symbols():
                for group in map_occupations:
                    if symb in group:
                        new_symbols.append(group[0])
                        break
        else:
            new_symbols = len(structure) * ["H"]
        structure_copy.set_chemical_symbols(new_symbols)
    dataset = spglib.get_symmetry_dataset(
        ase_atoms_to_spglib_cell(structure_copy), symprec=symprec
    )
    n_unitcells = np.linalg.det(dataset.transformation_matrix)

    equivalent_atoms = list(dataset.equivalent_atoms)
    wyckoffs = {}
    for index in set(equivalent_atoms):
        multiplicity = list(dataset.equivalent_atoms).count(index) / n_unitcells
        multiplicity = int(round(multiplicity))
        wyckoff = f"{multiplicity}{dataset.wyckoffs[index]}"
        if include_representative_atom_index:
            wyckoff += f"-{index}"
        wyckoffs[index] = wyckoff

    return [wyckoffs[equivalent_atoms[a.index]] for a in structure_copy]


"""
This module provides a simple wrapper of the ASE NeighborList class,
returning a list of Lattice Sites.
"""
from typing import Union
from ase import Atoms
from ase.neighborlist import NeighborList as ASENeighborList


def get_neighbor_lists(
    structure: Atoms, cutoffs: list, position_tolerance: float = 1e-5
) -> list[list[list[tuple]]]:
    """
    Returns a list of ASE neighbor lists given a configuration and cutoffs.

    Parameters
    ----------
    structure : Atoms
        ASE Atoms object representing the atomic configuration.
    cutoffs : list
        Positive floats indicating the cutoffs for different clusters.
    position_tolerance : float
        Tolerance applied when comparing positions in Cartesian coordinates.

    Returns
    -------
    neighbor_lists : list of list of list of tuples
        Outer list: one entry per cutoff.
        Middle list: one entry per atom.
        Inner list: neighbors of the atom as tuples (index, offset),
                    where 'offset' is the periodic cell offset vector.
    """
    if not isinstance(structure, Atoms):
        raise TypeError(f"Expected ASE Atoms object, got {type(structure)}")

    neighbor_lists = []

    for cutoff in cutoffs:
        nl = ASENeighborList(
            [cutoff / 2.0] * len(structure),
            skin=2 * position_tolerance,
            bothways=True,
            self_interaction=False,
        )
        nl.update(structure)

        cutoff_neighbors = []
        for i in range(len(structure)):
            indices, offsets = nl.get_neighbors(i)
            neighbors = [(idx, tuple(off)) for idx, off in zip(indices, offsets)]
            # Sort by offset first, then by index
            neighbors = sorted(neighbors, key=lambda x: (x[1], x[0]))
            cutoff_neighbors.append(neighbors)

        neighbor_lists.append(cutoff_neighbors)

    return neighbor_lists


import numpy as np
import spglib
from ase import Atoms


class MatrixOfEquivalentPositions:
    """
    This class handles a matrix of equivalent positions given the symmetry
    elements of an atomic structure.

    Note
    ----
    As a user you will usually not interact directly with objects of this type.

    Parameters
    ----------
    translations
        Translational symmetry elements.
    rotations
        Rotational symmetry elements.
    """

    def __init__(self, translations: np.ndarray, rotations: np.ndarray):
        if len(translations) != len(rotations):
            raise ValueError(
                f"The number of translations ({len(translations)})"
                f" must equal the number of rotations ({len(rotations)})."
            )
        self.n_symmetries = len(rotations)
        self.translations = np.array(translations)
        self.rotations = np.array(rotations)

    def build(self, fractional_positions: np.ndarray) -> None:
        """
        Builds a matrix of symmetry equivalent positions given a set of input
        coordinates using the rotational and translational symmetries provided upon
        initialization of the object.

        Parameters
        ----------
        fractional_positions
            Atomic positions in fractional coordinates.
            Dimensions: (number of atoms, 3 fractional coordinates).
        """
        positions = np.dot(self.rotations, fractional_positions.transpose())
        positions = np.moveaxis(positions, 2, 0)
        translations = self.translations[np.newaxis, :].repeat(
            len(fractional_positions), axis=0
        )
        positions += translations
        self.positions = positions

    def get_equivalent_positions(self) -> np.ndarray:
        """
        Returns the matrix of equivalent positions. Each row corresponds
        to a set of symmetry equivalent positions. The entry in the
        first column is commonly treated as the representative position.
        Dimensions: (number of atoms, number of symmetries, 3 fractional coordinates)
        """
        return self.positions


import numpy as np
import spglib
from ase import Atoms


class MatrixOfEquivalentPositions:
    """
    This class handles a matrix of equivalent positions given the symmetry
    elements of an atomic structure.

    Note
    ----
    As a user you will usually not interact directly with objects of this type.

    Parameters
    ----------
    translations
        Translational symmetry elements.
    rotations
        Rotational symmetry elements.
    """

    def __init__(self, translations: np.ndarray, rotations: np.ndarray):
        if len(translations) != len(rotations):
            raise ValueError(
                f"The number of translations ({len(translations)})"
                f" must equal the number of rotations ({len(rotations)})."
            )
        self.n_symmetries = len(rotations)
        self.translations = np.array(translations)
        self.rotations = np.array(rotations)

    def build(self, fractional_positions: np.ndarray) -> None:
        """
        Builds a matrix of symmetry equivalent positions given a set of input
        coordinates using the rotational and translational symmetries provided upon
        initialization of the object.

        Parameters
        ----------
        fractional_positions
            Atomic positions in fractional coordinates.
            Dimensions: (number of atoms, 3 fractional coordinates).
        """
        positions = np.dot(self.rotations, fractional_positions.transpose())
        positions = np.moveaxis(positions, 2, 0)
        translations = self.translations[np.newaxis, :].repeat(
            len(fractional_positions), axis=0
        )
        positions += translations
        self.positions = positions

    def get_equivalent_positions(self) -> np.ndarray:
        """
        Returns the matrix of equivalent positions. Each row corresponds
        to a set of symmetry equivalent positions. The entry in the
        first column is commonly treated as the representative position.
        Dimensions: (number of atoms, number of symmetries, 3 fractional coordinates)
        """
        return self.positions


import numpy as np
from ase import Atoms
from typing import List

Vector = List[float]


def get_fractional_positions_from_neighbor_list(
    structure: Atoms, neighbor_list: list
) -> List[Vector]:
    """
    Returns the fractional positions of the lattice sites in structure
    from a neighbor list, including periodic images.

    Parameters
    ----------
    structure : Atoms
        Input atomic structure.
    neighbor_list : list[list[tuple]]
        Each element is a list of neighbors, where each neighbor is
        a tuple (atom_index, image_vector).

    Returns
    -------
    List[Vector]
        Fractional positions of central atoms and their neighbors.
    """
    neighbor_positions = []

    for i, neighbors in enumerate(neighbor_list):
        # Add central atom
        neighbor_positions.append(structure.positions[i])

        # Add neighbor atoms with periodic image offsets
        for idx, image_vec in neighbors:
            image_vec = np.array(image_vec)  # e.g., [1, 0, -1]
            pos = structure.positions[idx] + np.dot(image_vec, structure.get_cell())
            neighbor_positions.append(pos)

    if len(neighbor_positions) > 0:
        # Convert to fractional coordinates
        fractional_positions = np.linalg.solve(
            structure.get_cell().T, np.array(neighbor_positions).T
        ).T
    else:
        fractional_positions = []

    return fractional_positions


def matrix_of_equivalent_positions_from_structure(
    structure: Atoms,
    cutoff: float,
    position_tolerance: float,
    symprec: float,
    find_primitive: bool = True,
):
    """Sets up a matrix of equivalent positions from an :class:`Atoms <ase.Atoms>` object.

    Parameters
    ----------
    structure
        Input structure.
    cutoff
        Cutoff radius.
    find_primitive
        If ``True`` the symmetries of the primitive structure will be employed.
    symprec
        Tolerance imposed when analyzing the symmetry using spglib.
    position_tolerance
        Tolerance applied when comparing positions in Cartesian coordinates.

    Returns
    -------
        The tuple that is returned comprises the matrix of equivalent positions,
        the primitive structure, and the neighbor list.
    """

    structure = structure.copy()
    structure_prim = structure
    if find_primitive:
        structure_prim = get_primitive_structure(structure, symprec=symprec)

    # get symmetry information
    structure_as_tuple = ase_atoms_to_spglib_cell(structure_prim)
    symmetry = spglib.get_symmetry(structure_as_tuple, symprec=symprec)
    translations = symmetry["translations"]
    rotations = symmetry["rotations"]

    # set up a MatrixOfEquivalentPositions object
    matrix_of_equivalent_positions = MatrixOfEquivalentPositions(
        translations, rotations
    )

    # create neighbor lists
    # prim_icet_structure = Structure.from_atoms(structure_prim)

    neighbor_list = get_neighbor_lists(
        structure_prim, [cutoff], position_tolerance=position_tolerance
    )[0]

    # print(neighbor_list)

    # get fractional positions for neighbor_list
    frac_positions = get_fractional_positions_from_neighbor_list(
        structure_prim, neighbor_list
    )

    if frac_positions is not None:
        matrix_of_equivalent_positions.build(frac_positions)

    return matrix_of_equivalent_positions, structure_prim, neighbor_list


import numpy as np
from ase import Atoms


def fractional_to_cartesian(
    fractional_coordinates: np.ndarray, cell: np.ndarray
) -> np.ndarray:
    """
    Convert fractional coordinates to Cartesian coordinates.

    Parameters
    ----------
    fractional_coordinates : np.ndarray, shape=(N,3)
        Fractional coordinates of atoms.
    cell : np.ndarray, shape=(3,3)
        Lattice vectors as rows.

    Returns
    -------
    cartesian_coordinates : np.ndarray, shape=(N,3)
        Cartesian coordinates.
    """
    return np.dot(fractional_coordinates, cell)


def find_atom_indices_by_positions(
    atoms: Atoms, positions: np.ndarray, tol: float = 1e-5
) -> list[int]:
    """
    Find atom indices in the ASE Atoms object that match given positions
    within a tolerance.

    Parameters
    ----------
    atoms : ase.Atoms
        ASE atoms object.
    positions : np.ndarray, shape=(N,3)
        Cartesian positions to match.
    tol : float
        Tolerance in Angstrom.

    Returns
    -------
    indices : list of int
        List of atom indices matching the positions.
    """
    indices = []
    for pos in positions:
        # Compute distance to all atoms
        distances = np.linalg.norm(atoms.positions - pos, axis=1)
        idx = np.where(distances < tol)[0]
        if len(idx) > 0:
            indices.append(int(idx[0]))
        else:
            raise ValueError(f"No atom found near position {pos}")
    return indices


def _fractional_to_cartesian(
    fractional_coordinates: list[list[float]], cell: np.ndarray
) -> list[float]:
    """
    Converts cell metrics from fractional to Cartesian coordinates.

    Parameters
    ----------
    fractional_coordinates
        List of fractional coordinates.
    cell
        Cell metric.
    """
    cartesian_coordinates = [np.dot(frac, cell) for frac in fractional_coordinates]
    return cartesian_coordinates


import numpy as np
from ase import Atoms
from typing import List, Tuple, Sequence


def find_lattice_sites_by_positions(
    structure: Atoms, positions: List[Sequence], fractional_position_tolerance: float
) -> List[Tuple[int, Tuple[int, int, int]]]:
    """
    Returns lattice sites (atom index + periodic image) corresponding to positions.

    Parameters
    ----------
    structure : Atoms
        ASE Atoms object with periodic boundary conditions.
    positions : list[Sequence]
        List of positions in Cartesian coordinates.
    fractional_position_tolerance : float
        Tolerance for matching positions in fractional coordinates.

    Returns
    -------
    List[Tuple[int, Tuple[int, int, int]]]
        List of lattice sites as (atom_index, image_vector) tuples.
    """
    lattice_sites = []
    cell = structure.get_cell()
    inv_cell = np.linalg.inv(cell.T)  # to convert Cartesian -> fractional

    for pos in positions:
        # Convert to fractional coordinates
        frac_pos = np.dot(inv_cell, np.array(pos))

        # Wrap into [0,1) for comparison
        frac_pos_wrapped = frac_pos % 1.0

        # Find the closest atom
        found = False
        for i, atom_frac in enumerate(structure.get_scaled_positions()):
            delta = frac_pos_wrapped - (atom_frac % 1.0)
            delta -= np.round(delta)  # account for periodic images
            if np.linalg.norm(delta) < fractional_position_tolerance:
                image_vec = tuple(np.round(frac_pos - atom_frac).astype(int))
                lattice_sites.append((i, image_vec))
                found = True
                break
        if not found:
            raise ValueError(f"Position {pos} not matched to any lattice site.")

    return lattice_sites


def _get_lattice_site_matrix_of_equivalent_positions(
    atoms: Atoms,
    matrix_of_equivalent_positions: MatrixOfEquivalentPositions,
    fractional_position_tolerance: float = 1e-5,
) -> list[list[int]]:
    """
    Transform a matrix of equivalent positions (fractional coordinates)
    into a matrix of lattice site indices.

    Parameters
    ----------
    atoms : ase.Atoms
        ASE atoms object representing the primitive cell.
    matrix_of_equivalent_positions : np.ndarray, shape=(M,N,3)
        Fractional positions for M rows of N atoms each.
    fractional_position_tolerance : float
        Tolerance in fractional coordinates.

    Returns
    -------
    lattice_site_matrix : list of list of int
        Matrix with the corresponding atom indices for each position.
    """

    eqpos_frac = matrix_of_equivalent_positions.get_equivalent_positions()
    eqpos_lattice_sites = []

    for row in eqpos_frac:
        positions = _fractional_to_cartesian(row, atoms.cell)
        lattice_sites = []
        if np.all(atoms.pbc):
            lattice_sites = find_lattice_sites_by_positions(
                atoms,
                positions=positions,
                fractional_position_tolerance=fractional_position_tolerance,
            )
        else:
            raise ValueError("Input structure must have periodic boundary conditions.")
        if lattice_sites is not None:
            eqpos_lattice_sites.append(lattice_sites)
        else:
            print(
                "Unable to transform any element in a column of the"
                " fractional matrix of equivalent positions to lattice site"
            )

    return eqpos_lattice_sites


from collections.abc import Iterable
import copy


def get_chemical_symbols(input_chemical_symbols, input_structure):
    """Returns chemical symbols using input structure and
    chemical symbols. Carries out multiple sanity checks."""

    # setup chemical symbols as List[List[str]]
    if all(isinstance(i, str) for i in input_chemical_symbols):
        chemical_symbols = [input_chemical_symbols] * len(input_structure)
    # also accept tuples and other iterables but not, e.g., List[List, str]
    # (need to check for str explicitly here because str is an Iterable)
    elif not all(
        isinstance(i, Iterable) and not isinstance(i, str)
        for i in input_chemical_symbols
    ):
        raise TypeError(
            "chemical_symbols must be List[str] or List[List[str]], not {}".format(
                type(input_chemical_symbols)
            )
        )
    elif len(input_chemical_symbols) != len(input_structure):
        msg = "chemical_symbols must have same length as structure. "
        msg += "len(chemical_symbols) = {}, len(structure)= {}".format(
            len(input_chemical_symbols), len(input_structure)
        )
        raise ValueError(msg)
    else:
        chemical_symbols = copy.deepcopy(input_chemical_symbols)

    for i, symbols in enumerate(chemical_symbols):
        if len(symbols) != len(set(symbols)):
            raise ValueError(
                "Found duplicates of allowed chemical symbols on site {}."
                " allowed species on  site {}= {}".format(i, i, symbols)
            )

    if len([tuple(sorted(s)) for s in chemical_symbols if len(s) > 1]) == 0:
        raise ValueError("No active sites found")

    return chemical_symbols


def GetInputsForOrbitList(structure, chemicalSymbols, cutoffs, symprec=1e-5):
    updateChemicalSymbols = get_chemical_symbols(chemicalSymbols, structure)

    occupied_primitive, primitive_chemical_symbols = get_occupied_primitive_structure(
        structure, updateChemicalSymbols, symprec=symprec
    )

    position_tolerance = symprec

    effective_box_size = abs(np.linalg.det(occupied_primitive.cell)) ** (1 / 3)
    tol = position_tolerance / effective_box_size
    tol = min(tol, position_tolerance / 5)
    fractional_position_tolerance = round(tol, -int(floor(log10(abs(tol)))))

    maxCutoff = np.max(cutoffs)

    matrix_of_equivalent_positions, prim_structure, _ = (
        matrix_of_equivalent_positions_from_structure(
            occupied_primitive,
            cutoff=maxCutoff,
            position_tolerance=position_tolerance,
            find_primitive=False,
            symprec=symprec,
        )
    )

    neighbor_lists = get_neighbor_lists(
        structure=prim_structure, cutoffs=cutoffs, position_tolerance=position_tolerance
    )

    pm_lattice_sites = _get_lattice_site_matrix_of_equivalent_positions(
        atoms=prim_structure,
        matrix_of_equivalent_positions=matrix_of_equivalent_positions,
        fractional_position_tolerance=fractional_position_tolerance,
    )
    # print(prim_structure)
    # print(pm_lattice_sites)
    # print(neighbor_lists)
    # print(position_tolerance)

    return (
        prim_structure,
        pm_lattice_sites,
        neighbor_lists,
        position_tolerance,
        fractional_position_tolerance,
    )


# ase atoms , list[list[int]], list[list[list[tuple]]], double


def GetInputsForOrbitListPrimStructre(
    structure, chemicalSymbols, cutoffs, symprec=1e-5
):
    prim_structure, _, _, _, _ = GetInputsForOrbitList(
        structure, chemicalSymbols, cutoffs, symprec
    )

    return prim_structure


def GetInputsForOrbitListPMLatticeSites(
    structure, chemicalSymbols, cutoffs, symprec=1e-5
):
    _, pm_lattice_sites, _, _, _ = GetInputsForOrbitList(
        structure, chemicalSymbols, cutoffs, symprec
    )

    return pm_lattice_sites


def GetInputsForOrbitListNeighbourList(
    structure, chemicalSymbols, cutoffs, symprec=1e-5
):
    _, _, neighbor_lists, _, _ = GetInputsForOrbitList(
        structure, chemicalSymbols, cutoffs, symprec
    )

    return neighbor_lists


def GetInputsForOrbitListTolerance(structure, chemicalSymbols, cutoffs, symprec=1e-5):
    _, _, _, position_tolerance, _ = GetInputsForOrbitList(
        structure, chemicalSymbols, cutoffs, symprec
    )

    return position_tolerance


def GetInputsForOrbitListFractionalTolerance(
    structure, chemicalSymbols, cutoffs, symprec=1e-5
):
    _, _, _, _, fractional_position_tolerance = GetInputsForOrbitList(
        structure, chemicalSymbols, cutoffs, symprec
    )

    return fractional_position_tolerance

