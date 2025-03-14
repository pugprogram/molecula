import networkx as nx
import matplotlib.pyplot as plt

def read_mol_file(filename):
    """
    Reads an OEChem mol file in V2000 format and returns the atom and bond
    information as separate lists.

    Args:
        filename (str): the name of the mol file

    Returns:
        atoms (list): a list of dictionaries representing the atoms in the molecule,
            with keys for element symbol, charge, and coordinates
        bonds (list): a list of tuples representing the bonds in the molecule,
            where each tuple contains the indices of the two atoms involved in the bond
            and the bond order
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[3].split()[0])
    num_bonds = int(lines[3].split()[1])
    atoms = []
    bonds = []
    for i in range(num_atoms):
        atom_line = lines[i+4].split()
        coords = tuple(map(float, [atom_line[0], atom_line[1], atom_line[2]]))
        atom = {'element': atom_line[3], 'charge': int(atom_line[4]), 'coords': coords}
        atoms.append(atom)
    for i in range(num_bonds):
        bond_line = lines[i+num_atoms+4].split()
        atom1 = int(bond_line[0])
        atom2 = int(bond_line[1])
        bond_order = int(bond_line[2])
        bonds.append((atom1, atom2, bond_order))
    return atoms, bonds
  
def generate_molecular_graph(atoms, bonds):
    """
    Generates a molecular graph from a list of atoms and bonds.

    Args:
        atoms (list): a list of dictionaries representing the atoms in the molecule,
            with keys for element symbol, charge, and coordinates
        bonds (list): a list of tuples representing the bonds in the molecule,
            where each tuple contains the indices of the two atoms involved in the bond
            and the bond order

    Returns:
        G (networkx.Graph): a NetworkX graph object representing the molecule,
            where each node represents an atom and each edge represents a bond,
            and each node and edge has attributes for element symbol and bond order, 
            respectively
    """
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        # add atom node to graph
        #Исключим водород сразу, чтобы потом не мучаться с ним
        if atom['element'] == 'H':
            continue
        node_label = str(i + 1)
        G.add_node(node_label, element=atom['element'])
        
    for bond in bonds:
        # add bond edge to graph
        atom1 = str(bond[0])
        atom2 = str(bond[1])
        #Исключим водород сразу, чтобы потом не мучаться с ним
        if atom1 in G.nodes and atom2 in G.nodes:
            bond_type = bond[2]
            G.add_edge(atom1, atom2, order=bond_type)
    return G
  
def canonical_Morgan(G):
    """
    Computes the canonical invariants and labels for each atom in a molecular graph.

    Args:
        G (networkx.Graph): a NetworkX graph object representing the molecule

    Returns:
        inv (dict): a dictionary of canonical invariants for each atom
        lab (dict): a dictionary of canonical labels for each atom
    """
    # initialize invariants and labels
    inv = {i: 1 for i in G.nodes()}
    lab = {i: 0 for i in G.nodes()}
    # compute invariants
    num_iterations = 1000
    flag = False
    inv, flag = compute_invariant_first_iteration(G, inv)
    if flag == False :
        
        for i in range (0, num_iterations):
            inv, flag = compute_invariant(G, inv)
            
            if flag == True:
                break

    if flag == False:
        print("Desicion not found!")
    else: 
        print("ALL OK")

    print(inv)

    lab = {node: G.nodes[node]['element'] for node in G.nodes()}

    print(lab)
    return inv, lab
  
def find_min_vertex(inv, lab):
    min_value = inv['1']
    min_item = '1'

    for key, value in inv.items():
        
        if value < min_value:
            min_item = key
            min_value = value
        elif value == min_value:
            current_priority = element_priority.get(lab[key], 0)
            min_priority = element_priority.get(lab[min_item], 0)
            
            if current_priority < min_priority:
                min_item = key
                min_value = value 


    return min_item

def compute_invariant(G, inv):
    inv_length = len(set(inv.values()))
    # compute invariant of each atom
    new_inv = {}

    for x in G.nodes():
        sum_vertex = 0

        for neighbor in G.neighbors(x):
            sum_vertex += inv[neighbor]
        
        new_inv[x] = inv[x] + sum_vertex

    new_inv_length = len(set(new_inv.values()))
    return new_inv, new_inv_length == inv_length

# Делаем отдельную функцию для 1 итерации, так как первая итерация не предполагает складывание самого узла графа
def compute_invariant_first_iteration(G, inv):
    inv_length = len(set(inv.values()))
    # compute invariant of each atom
    new_inv = {i: 0 for i in G.nodes()}
    
    for x in G.edges():
        new_inv[x[0]] += inv[x[1]]
        new_inv[x[1]] += inv[x[0]]

    new_inv_length = len(set(new_inv.values()))
    return new_inv, new_inv_length == inv_length
        
def find_lab (lab, inv, G, first_vertex):
    queue = []
    visited = set()
    lab[first_vertex] = 1
    queue.append(first_vertex)
    current_label = 1

    while queue:
        node = queue.pop(0)
        if node in visited:
            continue

        visited.add(node)
        lab[node] = current_label
        current_label += 1

        neighbors = sorted(
            G.neighbors(node),
            key = lambda x: (
                -inv[x]
            )
        )

        for neighbor in neighbors:
            if neighbor not in visited:
                queue.append(neighbor)

    

    return lab


value_order = {
    1 : "",
    2 : "=",
    3 : "#",
}

element_priority = {
    'C': 1,
    'N': 2,
    'O': 3,
    'S': 4, 
}

dict_for_cycle = {}
num_cycle = 0
neighbors_dict = {}

def find_smile_2(G, inv, lab, vertex, smile, visited):
    global neighbors_dict, num_cycle, dict_for_cycle, value_order, element_priority
    
    visited.add(vertex)  
    smile += lab[vertex]  # Добавляем метку вершины (атом) в строку

    neighbors = list(G.neighbors(vertex))  
    if not neighbors:
        return smile 

    # Сортируем соседей по приоритету обхода
    neighbors = sorted(
        neighbors,
        key=lambda x: (
            x not in visited,  # Непосещенные идут первыми
            -inv[x],  # Порядок обхода
            element_priority.get(G.nodes[x]['element'], 0)  # Приоритет элемента
        ),
        reverse=True
    )

    print(vertex)
    print(smile)
    for neighbor in neighbors:
        if (vertex == '4'):
            print("AAAA")

        if vertex == '5':
            print("SSS")
    
        if (neighbor in visited) and (neighbor in dict_for_cycle):
            smile += dict_for_cycle[neighbor]


        if (neighbor in visited) and (vertex in neighbors_dict[neighbor]):  
            num_cycle += 1  
            smile += str(num_cycle)  
            dict_for_cycle[neighbor] = num_cycle  
            continue  

        if (neighbor in visited):
            continue
        
        elif len(neighbors_dict[vertex]) == 1:  # Если нет ветвления
            connection = G.edges[vertex, neighbor]['order']
            smile += value_order[connection]  
            neighbors_dict[neighbor].remove(vertex)  
            neighbors_dict[vertex].remove(neighbor)  
            smile += find_smile(G, inv, lab, neighbor, "", visited)
            continue

        if (vertex == '5'):
            print("KEEE")
        neighbors_dict[neighbor].remove(vertex)  
        neighbors_dict[vertex].remove(neighbor)  
        new_str = find_smile(G, inv, lab, neighbor, "", visited)  

        smile += "("  
        connection = G.edges[vertex, neighbor]['order']  
        smile += value_order[connection]
        smile += new_str  
        smile += ")"  

    return smile 

beggin_vertex_symbol = {}
end_vertex_symbol = {}
#Смотрим на связи
relation = set()
num = 0

def find_smile(G, inv, lab, vertex, visited, smile, array_str, prev_vertex):
    global num, beggin_vertex_symbol, end_vertex_symbol
    visited.add(vertex)
    smile.append(vertex)

    smile_str = ""
    if prev_vertex != None:
        connection = G.edges[vertex, prev_vertex]['order']
        smile_str += value_order[connection]  
    smile_str += lab[vertex]
    array_str.append(smile_str)

    neighbors = list(G.neighbors(vertex))  
    if not neighbors:
        return smile, array_str 
    
    neighbors = sorted(
        neighbors,
        key=lambda x: (
            x not in visited,  # Непосещенные идут первыми
            -inv[x],  # Порядок обхода
            element_priority.get(G.nodes[x]['element'], 0)  # Приоритет элемента
        ),
        reverse=True
    )

    len_already_visited = 0
    for neighbor in neighbors:
        if neighbor in visited:
            len_already_visited += 1
    

    for neighbor in neighbors:
        #Если у нас есть связь с тем узлом, в котором мы уже были, то это явно указывает на цикл
        if (neighbor in visited) and ((str(neighbor) + str(vertex)) not in relation):
            num += 1
            beggin_vertex_symbol[neighbor] = ""
            end_vertex_symbol[neighbor] = ""
            beggin_vertex_symbol[vertex] = ""
            end_vertex_symbol[vertex] = str(num)
            relation.add(str(neighbor) + str(vertex))
            relation.add(str(vertex) + str(neighbor))


            for new in list(G.neighbors(neighbor)):
                if new == vertex: 
                    continue
                if new in beggin_vertex_symbol and (beggin_vertex_symbol[new] == '(' or beggin_vertex_symbol[new]>='1' and beggin_vertex_symbol[new]<=str(num)):
                    beggin_vertex_symbol[new] = ""
                    str_n = ""
                    str_n += str(num)
                    if new in end_vertex_symbol:
                        str_n += end_vertex_symbol[new]                    
                    end_vertex_symbol[new] = str_n
                    break

            continue

        # Тут рассмотрим случае ветвлений
        if (len(neighbors) - len_already_visited) > 1 :
            #Допускаем что у нас скобки всегда
            len_already_visited += 1
            beggin_vertex_symbol[neighbor] = '('
            #Добавляем связь по которой ходили
            relation.add(str(neighbor) + str(vertex))
            relation.add(str(vertex) + str(neighbor))
            smile, array_str = find_smile(G, inv, lab, neighbor, visited, smile, array_str, vertex)
            if beggin_vertex_symbol[neighbor] == '(':
                end_vertex_symbol[smile[-1]] = ')'
            continue
        
        if ((str(neighbor) + str(vertex)) not in relation):
        #главная ветка
            relation.add(str(neighbor) + str(vertex))
            relation.add(str(vertex) + str(neighbor))

            smile, array_str = find_smile(G, inv, lab, neighbor, visited, smile, array_str, vertex)

    return smile, array_str    

        

        
def find_smile_str(smile_array, smile_array_str):
    global beggin_vertex_symbol, end_vertex_symbol
    smile = ""
    

    for i in range(0, len(smile_array)):
        if smile_array[i] in beggin_vertex_symbol:
            smile += beggin_vertex_symbol[smile_array[i]]

        smile += smile_array_str[i]

        if smile_array[i] in end_vertex_symbol:
            smile += end_vertex_symbol[smile_array[i]]

    return smile

        






def create_dict_edges (G):
    global neighbors_dict
    neighbors_dict = {node: set() for node in G.nodes()}
    for u, v in G.edges():
        neighbors_dict[u].add(v)
        neighbors_dict[v].add(u)
    
    return 


def gen_smiles(molfile):
    """
    Generates a SMILES string from a mol file.

    Args:
        molfile (str): the name of the mol file

    Returns:
        smiles (str): the SMILES string
    """

    # Read mol file
    atoms, bonds = read_mol_file(molfile)

    # Generate molecular graph
    G = generate_molecular_graph(atoms, bonds)
    print(G)
    print("Узлы графа:", G.nodes())
    print("Ребра графа: ", G.edges())
    # Рисуем граф
    pos = nx.spring_layout(G)  # Позиционируем узлы с помощью spring layout
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=2000, font_size=16, font_weight='bold')

    # Добавляем метки узлов (элементы)
    labels = {node: G.nodes[node]['element'] for node in G.nodes()}
    print(labels)
    nx.draw_networkx_labels(G, pos, labels, font_size=16, font_color='black')

    # Добавляем метки ребер (порядок связи)
    edge_labels = {(u, v): G.edges[u, v]['order'] for u, v in G.edges()}
    print(edge_labels)
    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=16, font_color='red')

    # Отображаем граф
    plt.show()
    # Compute canonical invariants and labels

    create_dict_edges(G)
    print(neighbors_dict)
    inv, lab = canonical_Morgan(G)

    min_vertex = find_min_vertex(inv,lab)
    smile = []
    visited = set()
    array_str = []
    smile, array_str = find_smile(G, inv, lab, min_vertex, visited, smile, array_str, None)
    print("\n")
    print(beggin_vertex_symbol)
    print("\n")
    print(end_vertex_symbol)
    print("\n")

    print(array_str)


    smile_str = find_smile_str(smile, array_str)
    return smile_str


if __name__ == '__main__':
    molfile = '3-amino-2-naphthoic_acid.mol'
    smiles = gen_smiles(molfile)
    print(smiles)