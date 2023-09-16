#===================================================#
#                                                   #
# Program:  To visualize the structure from kwant   #
# Version:  1.0                                     #
# Author:   Wun-Hao Kang                            #
# Function: make_system                             #
# Class:    Bubel_Writer                            #
#                                                   #
#===================================================#

from numpy import array, zeros, linspace, logical_and, logical_or
from numpy import sqrt, sin, cos, arctan2, absolute, argmin
from math import pi, log10
from time import time

# Problem for lead: SHAPE_RANGE_3D does not support

def make_system(length: float, width: float, lead_width: float, staggered_potential: float, rashba: float):

    # ______
    # |     |_
    # |     :_ lead 0
    # |__.._|
    #    ||
    #  lead 1

    from kwant.builder import Builder, HoppingKind
    from kwant.lattice import Polyatomic, TranslationalSymmetry
    from numpy import array

    sigma_0 = array([[1, 0], [0, 1]], dtype=complex)
    sigma_x = array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = array([[0, -1j], [1j, 0]], dtype=complex)

    def scattering_region(pos):
        return -length/2 <= pos[0] <= length/2 and -width/2 <= pos[1] <= width/2
    def lead_0_region(pos):
        return -lead_width/2 <= pos[1] <= lead_width/2
    def lead_1_region(pos):
        return -lead_width/2 <= pos[0] <= lead_width/2
    def scattering_onsite(site, U):
        electric_term = -U*sigma_0 if site.pos[0]-site.pos[1] < 0 else U*sigma_0
        soc_term = (-1)**(int(site.family.name)%2)*staggered_potential*sigma_0
        return electric_term+soc_term
    def lead_0_onsite(site, U):
        electric_term = -U*sigma_0
        soc_term = (-1)**(int(site.family.name)%2)*staggered_potential*sigma_0
        return electric_term+soc_term
    def lead_1_onsite(site, U):
        electric_term = U*sigma_0
        soc_term = (-1)**(int(site.family.name)%2)*staggered_potential*sigma_0
        return electric_term+soc_term
    def hopping_energy(target, source):
        theta = arctan2(target.pos[1]-source.pos[1], target.pos[0]-source.pos[0])
        return 1j*2/3*rashba*(sigma_x*sin(theta)-sigma_y*cos(theta))

    start_time = time()

    lattice_spacing = 0.142 # unit: nm
    lattice_vectors = array([[3, 0], [0, sqrt(3)]], dtype=float)*lattice_spacing
    basis = array([[-1, 0], [-0.5, sqrt(3)/2], [0.5, sqrt(3)/2], [1, 0]], dtype=float)*lattice_spacing
    unit_cell = Polyatomic(lattice_vectors, basis, name=tuple(str(i) for i in range(basis.shape[0])), norbs=2)
    sublattices = unit_cell.sublattices
    hopping_relations = (((0, 0), 1, 0), ((0, 0), 2, 1), ((0, 0), 3, 2), ((1, 0), 0, 3), ((0, 1), 0, 1), ((0, 1), 3, 2))

    syst = Builder()
    syst[unit_cell.shape(scattering_region, (0, 0))] = scattering_onsite
    syst[[HoppingKind(hop[0], sublattices[hop[1]], sublattices[hop[2]]) for hop in hopping_relations]] = hopping_energy

    lead_0 = Builder(TranslationalSymmetry(unit_cell.vec((1, 0))))
    lead_0[unit_cell.shape(lead_0_region, (0, 0))] = lead_0_onsite
    lead_0[[HoppingKind(hop[0], sublattices[hop[1]], sublattices[hop[2]]) for hop in hopping_relations]] = hopping_energy
    syst.attach_lead(lead_0)

    lead_1 = Builder(TranslationalSymmetry(unit_cell.vec((0, -1))))
    lead_1[unit_cell.shape(lead_1_region, (0, 0))] = lead_1_onsite
    lead_1[[HoppingKind(hop[0], sublattices[hop[1]], sublattices[hop[2]]) for hop in hopping_relations]] = hopping_energy
    syst.attach_lead(lead_1)

    print(f"It took {time()-start_time:.2f} sec on building system.")

    return syst

class Bubel_Writer:
    def __init__(self, final_syst, basis_name: list=[]):
        # if you concern the permutation of the basis, please input the basis_name you want otherwise the permutation will be random.

        self.site_position = [] # the positions of each site
        self.site_name = [] # flag_0, the index of the name of site, bubel_viewer doesn't support non-digital flag
        self.site_part = [] # flag_1, 0: scattering region, 1: padding, >=2: lead
        self.site_num_orbs = [] # number of orbitals of each site
        self.site_neighbors = [] # the neighbors of each site
        self.translational_vectors = [] # the translational vector for each lead
        self.lead_frame = [] # the frame to show leads
        self.lead_site = [] # flag_2, use binary number to denote
        self.from_atom_ids = [] # be used for hamiltonian
        self.to_atom_ids = [] # be used for hamiltonian
        self.from_orbital = [] # be used for hamiltonian
        self.to_orbital = [] # be used for hamiltonian
        self.__index_transform = [] # index of created atom -> index of site in each lead
        self.__read_info(final_syst, basis_name) # the information of lattice

    def __read_info(self, final_syst, basis_name):
        sites = final_syst.sites

        # scattering region
        self.site_position.extend([site.pos for site in sites])
        self.site_part.extend([0]*len(sites))
        self.site_num_orbs.extend([site.family.norbs for site in sites])
        temp_neighbors = list(map(final_syst.graph.out_neighbors, range(len(self.site_position))))
        self.site_neighbors.extend(list(map(list, temp_neighbors))) # the type of final_syst.graph.out_neighbors() is not list, it's Slice
        for site in sites:
            if site.family.name in basis_name:
                self.site_name.append(basis_name.index(site.family.name))
            else:
                self.site_name.append(len(basis_name))
                basis_name.append(site.family.name)
        del temp_neighbors
        self.lead_site.extend([0]*len(sites))

        # leads
        for i in range(len(final_syst.leads)):
            temp_sites = final_syst.leads[i].sites

            self.translational_vectors.append([final_syst.leads[i].symmetry.periods[0, j] for j in range(final_syst.leads[i].symmetry.periods.shape[1])])
            for index in final_syst.lead_paddings[i]:
                # the atoms of lead padding are exist in the scattering region but this is not what we defined.
                self.site_part[index] = 1

            # To created virtual atom for the lead
            temp_shift = [] # (n1, n2, n3) for shifting site in lead to the position we want. Assume the relation between the fist site in interface and the first site in the lead is integer times translational vector
            for j in range(len(self.translational_vectors[-1])):
                temp_shift.append(abs(sites[final_syst.lead_interfaces[i][0]].tag[j]-final_syst.leads[i].sites[0].tag[j]))
                temp_shift [-1] += 0 if temp_shift[-1] == 0 else 1

            temp_global_index = []; temp_frame = []
            temp_interface_pos = array([sites[site_i].pos for site_i in final_syst.lead_interfaces[i]])
            start_index = len(self.site_position); current_global_index = len(self.site_position);
            for j, site in enumerate(temp_sites):
                temp_new_pos = [site.pos[k]+temp_shift[k]*self.translational_vectors[-1][k] for k in range(len(temp_shift))]
                isinterface = sum([(temp_interface_pos[:, k]-temp_new_pos[k])**2 for k in range(len(temp_shift))]) < 10**(-16)
                if isinterface.sum() == 0:
                    temp_global_index.append(current_global_index)
                    self.site_position.append(temp_new_pos)
                    self.site_name.append(basis_name.index(temp_sites[j].family.name))
                    self.site_part.append(2+i)
                    self.site_num_orbs.append(temp_sites[j].family.norbs)
                    self.lead_site.append(0)
                    current_global_index += 1
                else:
                    temp_global_index.append(final_syst.lead_interfaces[i][isinterface][0])
                # lead frame
                if j == 0:
                    temp_frame.extend(temp_new_pos)
                    temp_frame.extend(temp_new_pos)
                else:
                    for k in range(len(temp_new_pos)):
                        temp_frame[k] = temp_new_pos[k] if temp_new_pos[k] < temp_frame[k] else temp_frame[k]
                        temp_frame[len(temp_new_pos)+k] = temp_new_pos[k] if temp_new_pos[k] > temp_frame[len(temp_new_pos)+k] else temp_frame[len(temp_new_pos)+k]
                # lead site for flag_2
                temp_site_pos = array(self.site_position, dtype=float)
                for k in range(len(temp_new_pos)):
                    temp_site_pos[:, k] -= site.pos[k]
                self.lead_site[argmin((temp_site_pos**2).sum(axis=1))] += 2**(i)
            self.lead_frame.append(temp_frame)

            for j in range(len(temp_sites)):
                temp_neighbors = list(final_syst.leads[i].graph.out_neighbors(j))
                if temp_global_index[j] < start_index:
                    # interface
                    for index in temp_neighbors:
                        if temp_global_index[index] >= start_index:
                            # add lead atom to original neighbors
                            self.site_neighbors[temp_global_index[j]].append(temp_global_index[index])
                    self.site_neighbors[temp_global_index[j]].sort()
                else:
                    # lead
                    for k, index in enumerate(temp_neighbors):
                        temp_neighbors[k] = temp_global_index[index]
                    self.site_neighbors.append(sorted(temp_neighbors))
            self.__index_transform.append(temp_global_index)

        del temp_new_pos, temp_frame, temp_interface_pos, temp_global_index, temp_neighbors

        # common part
        for i in range(len(self.site_position)):
            temp_nbs = sorted([i]+list(self.site_neighbors[i]))
            for j in range(len(temp_nbs)):
                self.from_atom_ids.extend([i]*self.site_num_orbs[i]*self.site_num_orbs[j])
                self.to_atom_ids.extend([temp_nbs[j]]*self.site_num_orbs[i]*self.site_num_orbs[j])
                self.from_orbital.extend([int(k/self.site_num_orbs[j]) for k in range(self.site_num_orbs[i]*self.site_num_orbs[j])])
                self.to_orbital.extend(list(range(self.site_num_orbs[j]))*self.site_num_orbs[i])
        del temp_nbs

        self.site_position = array(self.site_position, dtype=float)
        self.site_part = array(self.site_part, dtype=int)
        self.from_atom_ids = array(self.from_atom_ids, dtype=int)
        self.to_atom_ids = array(self.to_atom_ids, dtype=int)
        self.translational_vectors = array(self.translational_vectors, dtype=float)

    def __float2str(self, num: float):
        text = "" if num < 0 else " "
        text += f"{num:.8e}"
        return text
    def __int2str(self, num: float, num_digits: int):
        text = str(num)
        text = " "*(num_digits-len(text))+text
        return text
    def __get_hamiltonian(self, final_syst, *args, params: dict=None):
        h = zeros((self.from_atom_ids.size, 2), dtype=float)
        index = 0
        while index < self.from_atom_ids.size:
            isnolead = self.site_part[self.from_atom_ids[index]] < 2 and self.site_part[self.to_atom_ids[index]] < 2
            if isnolead:
                # part 0,1 <-> part 0,1
                temp_h = final_syst.hamiltonian(self.to_atom_ids[index], self.from_atom_ids[index], *args, params=params)
            else:
                # part 0,1 <-> part 2,3,... or part 2,3,.. <-> part 2,3,...
                lead_index = self.site_part[self.from_atom_ids[index]]-2 if self.site_part[self.from_atom_ids[index]] >= 2 else self.site_part[self.to_atom_ids[index]]-2
                new_from_atom_ids = self.__index_transform[lead_index].index(self.from_atom_ids[index])
                new_to_atom_ids = self.__index_transform[lead_index].index(self.to_atom_ids[index])
                temp_h = final_syst.leads[lead_index].hamiltonian(new_to_atom_ids, new_from_atom_ids, *args, params=params)
            # Though I think column is the source and row is the target and the matrix should be flatten in the order of "F", but the bubel viewer says I need to use the order of "C"
            h[index:index+temp_h.size, 0] = temp_h.flatten(order="C").real
            h[index:index+temp_h.size, 1] = temp_h.flatten(order="C").imag
            index += temp_h.size
        # To remvoe -0.00000000e+00
        mask = absolute(h[:, 0]) < 10**(-16); h[mask, 0] = 0;
        mask = absolute(h[:, 1]) < 10**(-16); h[mask, 1] = 0;
        return h

    def write(self, fname: str, final_syst, *args, params: dict=None):
        text = "<system>\n"
        text += "  <lattice>\n"
        text += "    <atoms>\n"

        # atom's information
        for i in range(self.site_position.shape[0]):
            text += "      <d>"
            # POS_X, POS_Y, POS_Z
            text += f" {self.__float2str(self.site_position[i, 0])} {self.__float2str(self.site_position[i, 1])}"
            if self.site_position.shape[1] < 3:
                text += f" {0.0:.8e}"
            else:
                text += f" {self.__float2str(self.site_position[i, 2])}"
            # FLAG_0
            text += f" {self.__int2str(self.site_name[i], int(log10(max(self.site_name)))+1)}"
            # FLAG_1
            text += f" {self.__int2str(self.site_part[i], int(log10(self.site_part.max()))+1)}"
            # FLAG_2
            text += f" {self.__int2str(self.lead_site[i], int(log10(max(self.lead_site)))+1)}"
            # ACTIVE: no idea
            text += " 1"
            # NUM_ORBITALS
            text += f" {self.__int2str(self.site_num_orbs[i], int(log10(max(self.site_num_orbs))+1))}"
            # NUM_BOUNDS
            text += f" {self.__int2str(len(self.site_neighbors[i]), int(log10(max(list(map(len, self.site_neighbors)))))+1)}"
            text += " </d>\n"
        text += "    </atoms>\n"

        # Connections: onsite and hoppings
        hamiltonian = self.__get_hamiltonian(final_syst, *args, params=params)
        text += "    <connections>\n"
        for i in range(self.from_atom_ids.size):
            # FROM_ATOM_ID
            text += f"      <d> {self.__int2str(self.from_atom_ids[i]+1, int(log10(self.from_atom_ids.max()+1))+1)}"
            # TO_ATOM_ID
            text += f" {self.__int2str(self.to_atom_ids[i]+1, int(log10(self.to_atom_ids.max()+1))+1)}"
            # FROM_ORBITAL_I
            text += f" {self.__int2str(self.from_orbital[i]+1, int(log10(max(self.from_orbital)+1))+1)}"
            # TO_ORBITAL_J
            text += f" {self.__int2str(self.to_orbital[i]+1, int(log10(max(self.to_orbital)+1))+1)}"
            # RE_VALUE IM_VALUE
            text += f" {self.__float2str(hamiltonian[i, 0])} {self.__float2str(hamiltonian[i, 1])} </d>\n"
        text += "    </connections>\n"
        text += "  </lattice>\n"

        # Lead
        for i in range(self.site_part.max()-1):
            text += "  <lead>\n"
            # Lead data
            text += "    <shape_type>SHAPE_RECTANGLE_XY</shape_type>\n"
            text += "    <shape_data>\n"
            # X_MIN Y_MIN X_MAX Y_MAX
            text += f"      {self.__float2str(self.lead_frame[i][0])} {self.__float2str(self.lead_frame[i][1])} {self.__float2str(self.lead_frame[i][2])} {self.__float2str(self.lead_frame[i][3])}\n"
            text += "    </shape_data>\n"
            # Vector
            text += "    <vector>\n"
            text += f"      {self.__float2str(self.translational_vectors[i][0])} {self.__float2str(self.translational_vectors[i][1])} {0.0:.8e}\n"
            text += "    </vector>\n"

            lead_atom_index = array(self.__index_transform[i], dtype=int)
            lead_atom_part = self.site_part[lead_atom_index]
            mask = lead_atom_part == i+2
            # Atoms: The interface atoms
            text += "    <atoms>\n"
            for j in lead_atom_index[~mask]:
                text += f"      <d> {self.__int2str(j+1, int(log10(self.site_position.shape[0]+1)+1))} </d>\n"
            text += "    </atoms>\n"
            # Next_atoms: the atoms in lead but exclude the interface atoms
            text += "    <next_atoms>\n"
            for j in lead_atom_index[mask]:
                text += f"      <d> {self.__int2str(j+1, int(log10(self.site_position.shape[0]+1)+1))} </d>\n"
            text += "    </next_atoms>\n"

            indices = linspace(0, self.from_atom_ids.size-1, self.from_atom_ids.size, dtype=int)
            isleadfrom = zeros(self.from_atom_ids.size, dtype=bool); isleadto = isleadfrom.copy();
            for j in lead_atom_index:
                isleadfrom = logical_or(isleadfrom, self.from_atom_ids == j)
                isleadto = logical_or(isleadto, self.to_atom_ids == j)
            isleadfrom = logical_and(isleadfrom, isleadto) # atoms in lead
            mask2 = zeros(self.from_atom_ids.size, dtype=bool); mask3 = mask2.copy();
            for j in lead_atom_index[~mask]:
                mask2 = logical_or(mask2, self.from_atom_ids == j)
                mask3 = logical_or(mask3, self.to_atom_ids == j)
            mask2 = logical_and(logical_or(mask2, mask3), ~logical_and(mask2, mask3)) # interface atoms in lead exculde hopping of interface atoms
            # Lead_coupling: intercell hopping
            text += "    <lead_coupling>\n"
            for j in indices[logical_and(mask2, isleadfrom)]:
                # FROM_ATOM_ID
                text += f"      <d> {self.__int2str(self.from_atom_ids[j]+1, int(log10(self.from_atom_ids.max()+1))+1)}"
                # TO_ATOM_ID
                text += f" {self.__int2str(self.to_atom_ids[j]+1, int(log10(self.to_atom_ids.max()+1))+1)}"
                # FROM_ORBITAL_I
                text += f" {self.__int2str(self.from_orbital[j]+1, int(log10(max(self.from_orbital)+1))+1)}"
                # TO_ORBITAL_J
                text += f" {self.__int2str(self.to_orbital[j]+1, int(log10(max(self.to_orbital)+1))+1)} </d>\n"
            text += "    </lead_coupling>\n"

            mask2 = zeros(self.from_atom_ids.size, dtype=bool); mask3 = mask2.copy();
            for j in lead_atom_index[mask]:
                mask2 = logical_or(mask2, self.from_atom_ids == j)
                mask3 = logical_or(mask3, self.to_atom_ids == j)
            mask2 = logical_and(mask2, mask3) # atoms in the lead exclude the interface atoms
            # Inner_coupling: intracell hopping
            text += "    <inner_coupling>\n"
            for j in indices[logical_and(mask2, isleadfrom)]:
                # FROM_ATOM_ID
                text += f"      <d> {self.__int2str(self.from_atom_ids[j]+1, int(log10(self.from_atom_ids.max()+1))+1)}"
                # TO_ATOM_ID
                text += f" {self.__int2str(self.to_atom_ids[j]+1, int(log10(self.to_atom_ids.max()+1))+1)}"
                # FROM_ORBITAL_I
                text += f" {self.__int2str(self.from_orbital[j]+1, int(log10(max(self.from_orbital)+1))+1)}"
                # TO_ORBITAL_J
                text += f" {self.__int2str(self.to_orbital[j]+1, int(log10(max(self.to_orbital)+1))+1)} </d>\n"
            text += "    </inner_coupling>\n"
            text += "  </lead>\n"
        text += "</system>\n"

        with open(f"{fname}.xml", "w") as f:
            f.write(text)

if __name__ == "__main__":

    from matplotlib import pyplot

    ispadding = True
    length_factor = 4 if ispadding else 5 # number of unit cell along x-axis.
    width_factor = 3 if ispadding else 4 # number of unit cell along y-axis.
    lead_factor = 2 # number of unit cell along the width of lead_0

    length = length_factor*3*0.142 # unit: nm
    width = width_factor*sqrt(3)*0.142 # unit: nm
    lead_width = lead_factor*3*0.142 # unit: nm

    staggered_potential = 0.001 # unit: eV
    rashba = 0.01 # unit: eV

    device = make_system(length, width, lead_width, staggered_potential, rashba)

    #============================================================================================================================================#
    #                                                           Structure                                                                        #
    #============================================================================================================================================#
    if True:
        from kwant.plotter import plot as kwant_plot

        fig, ax = pyplot.subplots(dpi=300)
        kwant_plot(device, ax=ax)
        ax.vlines((-length/2, length/2), -width/2, width/2, color="k", linewidth=1, linestyles="dashed")
        ax.hlines((-width/2, width/2), -length/2, length/2, color="k", linewidth=1, linestyles="dashed")

        ax.set_xlabel("x (nm)")
        ax.set_xticks((-length/2, 0, length/2))
        ax.set_ylabel("y (nm)")
        ax.set_yticks((-width/2, 0, width/2))
        ax.set_aspect("equal")

        pyplot.tight_layout()
        pyplot.savefig(f"./Graphene_ribbon_{'include_padding' if ispadding else 'no_padding'}.png")

    #============================================================================================================================================#
    #                                                             Output                                                                         #
    #============================================================================================================================================#
    if True:

        final_syst = device.finalized()
        U = 0.1

        kwant_syst = Bubel_Writer(final_syst, ("0", "1", "2", "3"))
        kwant_syst.write(f"bubel_viewer_graphene_ribbon_{'include_padding' if ispadding else 'no_padding'}", final_syst, params=dict(U=U))
