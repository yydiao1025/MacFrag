
# MacFrag is an efficient molecule fragmentation method, which is capable of segmenting large-scale molecules in a rapid speed and 
# generating diverse fragments that are more compliant with the ¡°Rule of Three¡±. 

# MacFrag is developed and maintained by Prof. HongLin Li's Group, School of Pharmacy, East China University of Science & Technology, Shanghai 200237, China. 
# http://lilab-ecust.cn/



from rdkit import Chem
from rdkit.Chem.BRICS import BreakBRICSBonds
from rdkit import RDLogger

import copy
from itertools import combinations
#from simple import simple_iter
from igraph import Graph
import argparse

RDLogger.DisableLog('rdApp.*')


def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-input_file', '-i', required=True,                        
                        help='.smi or .sdf file of molecules to be fragmented')    
    
    parser.add_argument('-output_path', '-o', required=True,
                        help='path of the output fragments file')
    
    parser.add_argument('-maxBlocks', required=True,
                        help='the maximum number of building blocks that the fragments contain')
    
    parser.add_argument('-maxSR', required=True,
                        help='only cyclic bonds in smallest SSSR ring of size larger than this value will be cleaved')

    parser.add_argument('-asMols', required=True,
                        help='True of False; if True, MacFrag will reture fragments as molecules and the fragments.sdf file will be output; if False, MacFrag will reture fragments.smi file with fragments representd as SMILES strings')
    
    parser.add_argument('-minFragAtoms', required=True,
                        help='the minimum number of atoms that the fragments contain')

    return parser.parse_args()


# SMARTS representations of the atomic environments for molecule fragmentation.
environs = {
  'L1': '[C;D3]([#0,#6,#7,#8])(=O)',  #original L1
  'L2': '[O;D2]-[#0,#6,#1]',  #original L3
  'L3': '[C;!D1;!$(C=*)]-[#6]',  #original L4
  'L4': '[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]',  #original L5
  'L5': '[C;D2,D3]-[#6]',  #original L7
  'L6': '[C;!D1;!$(C!-*)]',  #original L8
  'L61': '[C;R1;!D1;!$(C!-*)]',
  'L7': '[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]',  #original L9
  'L8': '[N;R;$(N(@C(=O))@[#6,#7,#8,#16])]',  #original L10
  'L9': '[S;D2](-[#0,#6])',  #original L11
  'L10': '[S;D4]([#6,#0])(=O)(=O)',  #original L12
  
  'L11': '[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]',  #original L13
  'L111': '[C;R2;$(C(-;@[C,N,O,S])-;@[N,O,S])]',
  'L112': '[C;R1;$(C(-;@[C,N,O,S;R2])-;@[N,O,S;R2])]',
  
  'L12': '[c;$(c(:[c,n,o,s]):[n,o,s])]',  #original L14
  
  'L13': '[C;$(C(-;@C)-;@C)]',  #original L15
  'L131': '[C;R2;$(C(-;@C)-;@C)]',
  'L132': '[C;R1;$(C(-;@[C;R2])-;@[C;R2])]',
  
  'L14': '[c;$(c(:c):c)]', #original L16
}
reactionDefs = (
  # L1
   [('1', '2', '-'),
    ('1', '4', '-'),
    ('1', '8', '-'),
    ('1', '11', '-'),
    ('1', '12', '-'),
    ('1', '13', '-'),
    ('1', '14', '-')],

  # L2
    [('2', '3', '-'),
    ('2', '11', '-'),
    ('2', '12', '-'),
    ('2', '13', '-'),
    ('2', '14', '-')],

  # L3
    [('3', '4', '-'),
    ('3', '9', '-')],

  # L4
    [('4', '10', '-'),
    ('4', '12', '-'),
    ('4', '14', '-'),
    ('4', '11', '-'),
    ('4', '13', '-')],

  # L5
    [('5', '5', '=')],

  # L6
    [('6', '7', '-'),
    ('6', '8', '-'),
    ('6', '11', '-;!@'),
    ('6', '12', '-'),
    ('6', '13', '-;!@'),
    ('6', '14', '-')],
     
  # L61
    [('61', '111', '-;@'),
    ('61', '131', '-;@')],   

  # L7
    [('7', '11', '-'),  
    ('7', '12', '-'),  
    ('7', '13', '-'),
    ('7', '14', '-')],

  # L8
    [('8', '11', '-'),
    ('8', '12', '-'),
    ('8', '13', '-'),
    ('8', '14', '-')],

  # L9
    [('9', '11', '-'),
    ('9', '12', '-'),
    ('9', '13', '-'),
    ('9', '14', '-')],

  # L11
    [('11', '12', '-'),
    ('11', '13', '-;!@'),
    ('11', '14', '-')],
     
  # L112
    [('112', '132', '-;@')],

  # L12
    [('12', '12', '-'),  
    ('12', '13', '-'),
    ('12', '14', '-')],

  # L13
    [('13', '14', '-')],

  # L14
    [('14', '14', '-')],  
    )

environMatchers = {}
for env, sma in environs.items():
    
    environMatchers[env] = Chem.MolFromSmarts(sma)
   
        
bondMatchers = []

for compats in reactionDefs:
    
    tmp = []
    for i1, i2, bType in compats:
        e1 = environs['L%s' % i1]
        e2 = environs['L%s' % i2]
        patt = '[$(%s)]%s[$(%s)]' % (e1, bType, e2)
        patt = Chem.MolFromSmarts(patt)
        tmp.append((i1, i2, bType, patt))
    bondMatchers.append(tmp)


def SSSRsize_filter(bond,maxSR=8):
    
    judge=True
    for i in range(3,maxSR+1):
        if bond.IsInRingSize(i) :
            judge=False
            break           
    return judge
    
    
def searchBonds(mol,maxSR=8):
    
    bondsDone = set()
    
    envMatches = {}   
    for env, patt in environMatchers.items():
        envMatches[env] = mol.HasSubstructMatch(patt)
        
    for compats in bondMatchers:
                
        for i1, i2, bType, patt in compats:
            if not envMatches['L' + i1] or not envMatches['L' + i2]:
                continue
            
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                if match not in bondsDone and (match[1], match[0]) not in bondsDone:
                    bond=mol.GetBondBetweenAtoms(match[0],match[1])
                    
                    if not bond.IsInRing():
                        bondsDone.add(match) 
                        yield (((match[0], match[1]), (i1, i2)))
                    elif bond.IsInRing() and  SSSRsize_filter(bond,maxSR=maxSR):
                        bondsDone.add(match)
                        yield (((match[0], match[1]), (i1, i2)))
                    
def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def mol_remove_atom_mapnumber(mol):
    for atom in mol.GetAtoms():
        atom.ClearProp('molAtomMapNumber')
    return mol

def Get_block_index(blocks):   
    
    block_index={}
    i=0
    for bs in blocks:
        tmp=[a.GetAtomMapNum() for a in bs.GetAtoms() if a.GetSymbol()!='*']
#        tmp=list(set(tmp))
        block_index[tuple(tmp)]=i
        i+=1      
    return block_index

def extrac_submol(mol,atomList,link):
    aList_mol=list(range(mol.GetNumAtoms()))
    aList_link=list(set([a[0] for a in link]))
    aList_submol=list(set(atomList+aList_link))
    aList_remove=[a for a in aList_mol if a not in aList_submol]
    eMol=Chem.RWMol(mol)
    
    aList_bbond=[a for a in aList_link if a not in atomList]
    for b in combinations(aList_bbond,2):
        eMol.RemoveBond(b[0],b[1])
    
    aList_remove.sort(reverse=True)
    for a in aList_remove:
        eMol.RemoveAtom(a) 

            
    for ba,btype in link:
        if ba in atomList:
            continue
        tmpatom=[a for a in eMol.GetAtoms() if a.GetAtomMapNum()==ba][0]
        tmpatom.SetIsAromatic(False)
        tmpatom.SetAtomicNum(0)
        tmpatom.SetIsotope(int(btype))
        tmpatom.SetNoImplicit(True)
            
    frag=eMol.GetMol()
   
    for a in frag.GetAtoms():
        a.ClearProp('molAtomMapNumber')   
    return frag
 
    
def simple_iter(graph, k):
    cis=[]
#    numb_at_most_k = 0
    # colors for vertices (igraph is to slow)
    colors = graph.vcount()*[-1]
    # consider only graphs with at least k vertices
    for i in range(graph.vcount() -1, -1, -1):
        # subgraph_set
        subgraph_set = [i]
        # print first induced subgraph conisting only of vertex i
#        subgraph_file.write(print_names(graph, subgraph_set))
        cis.append(copy.deepcopy(subgraph_set))
#        numb_at_most_k += 1
        # extension set
        extension = []
        # list of lists of all exclusive neighbors
        ex_neighs = [[]]
        # color closed neighborhood of start vertex
        colors[i] = 0
        for vertex in graph.neighbors(i):
            if vertex >= i:
                break
            colors[vertex] = 1
            extension.append(vertex)
            ex_neighs[0].append(vertex)
        # lists for the pointers of the extension set
        pointers = [extension.__len__()-1]
        # we save the index in whose closed neighborhood the actual pointers are, to perform the jumps
        poi_col = [1]
        # subgraph size
        sub_size = 1
        # enumerate all subgraphs containing vertex i
        while subgraph_set != []:
            # for each vertex in the extension set, create a new branch
            # if the actual pointer points to null (-1), then the corresponding extension set is empty
            while pointers[-1] > -1:
                last = pointers[-1]
                vertex = extension[last]
                ver_col = colors[vertex]
                # move pointer one to the left
                pointers[-1] -= 1
                act_vertex = pointers[-1]
                # check if vertex is an exclusive neighbor of the vertex last added to the subgraph set
                if ver_col == sub_size:
                    # if yes delete this vertex from the extension set
                    extension.pop()
                # check if pointer jumps, in other words check if the new next vertex in the extension set has a different color than the current vertex
                if act_vertex > -1:
                    if colors[extension[act_vertex]] < poi_col[-1]:
                        # jump with pointer (we stay at the same point if the color of the actual vertex is one less)
                        pointers[-1] = pointers[poi_col[-1]-2]
                        # update color
                        poi_col[-1] = colors[extension[pointers[-1]]]
                # create next child
                subgraph_set.append(vertex)
#                numb_at_most_k += 1
#                subgraph_file.write(print_names(graph, subgraph_set))
                cis.append(copy.deepcopy(subgraph_set))
                sub_size += 1
                # check if adding this vertex leads to a solution
                if sub_size == k:
                    # check time limit
#                    time_now = clock()
#                    if time_now > time_max:
 #                       return numb_at_most_k
                    subgraph_set.pop()
                    sub_size -= 1
                # otherwise create new enumeration tree node
                else:
                    # find the exclusive neighbors of vertex
                    ex_neighs.append([])
                    found = False
                    for neig in graph.neighbors(vertex):
                        if neig >= i:
                            break
                        if colors[neig] == -1:
                            colors[neig] = sub_size
                            ex_neighs[-1].append(neig)
                            extension.append(neig)
                            found = True
                    # if there is an exclusive neighbor, the new pointer points to him, otherwise to the old pointer
                    if found == True:
                        pointers.append(extension.__len__()-1)
                        poi_col.append(sub_size)
                    else:
                        pointers.append(pointers[-1])
                        poi_col.append(poi_col[-1])
            # now the last pointer points to null, so restore the parent
            pointers.pop()
            poi_col.pop()
            subgraph_set.pop()
            for vertex in ex_neighs[-1]:
                colors[vertex] = -1
            ex_neighs.pop()
            sub_size -= 1
        # remove color from the root
        colors[i] = -1
    return cis       

       
def MacFrag(mol,maxBlocks=6,maxSR=8,asMols=False,minFragAtoms=1):
    
    fragPool={}
    
    mol=mol_with_atom_index(mol)    
   
    bonds=list(searchBonds(mol,maxSR=maxSR))
    fragmol=BreakBRICSBonds(mol,bonds=bonds)
    blocks=Chem.GetMolFrags(fragmol,asMols=True)
    for block in blocks:
        tmp=copy.deepcopy(block)
        tmp_smiles=Chem.MolToSmiles(tmp)  
        nAtoms=tmp.GetNumAtoms(onlyExplicit=True)
        if nAtoms - tmp_smiles.count('*') >= minFragAtoms:
            for a in tmp.GetAtoms():
                a.ClearProp('molAtomMapNumber') 
            fragPool[tmp]=Chem.MolToSmiles(tmp)

        
    if maxBlocks>len(blocks):
        maxBlocks=len(blocks)
    if maxBlocks==1:
        if asMols:
            return fragPool.keys()
        else:
            return list(set(fragPool.values()))
       
    block_index=Get_block_index(blocks)
    index_block={block_index[key]:key for key in block_index.keys()}
                

    bond_block={}
    point=False
    block_link={block:set() for block in block_index.keys()}
    for b in bonds:
        ba1,ba2=b[0]
        btype1,btype2=b[1]
        for block in block_index.keys():
            if ba1 in block:
                bond_block[ba1]=block
                block_link[block].add((ba2,btype2))
                point=True
            if ba2 in block:
                bond_block[ba2]=block
                block_link[block].add((ba1,btype1))
                point=True
            if point==True:continue
    
    n=len(index_block.keys())    
    edges=[]
    for b in bonds:
        ba1,ba2=b[0]
        edges.append((block_index[bond_block[ba1]],block_index[bond_block[ba2]]))
    
    graph=Graph(n=n,edges=edges,directed=False)
#    all_cis=rssp(graph,k=maxBlocks)
    all_cis=simple_iter(graph,k=maxBlocks)
    for i in range(n):
        all_cis.remove([i])
    
    sub_link={}
    for cis in all_cis:
        tmp=[]
        tmp2=set()
        for ni in cis:
            tmp.extend(list(index_block[ni]))
            tmp2=tmp2.union(block_link[index_block[ni]])
        sub_link[tuple(tmp)]=tmp2  

  
    for fa in sub_link:
        frag=extrac_submol(mol,atomList=list(fa),link=sub_link[fa])
        frag_smiles=Chem.MolToSmiles(frag)    
        nAtoms=frag.GetNumAtoms(onlyExplicit=True)
        if nAtoms - frag_smiles.count('*') >= minFragAtoms:
            fragPool[frag]=frag_smiles
                
        
    if asMols:
        return list(fragPool.keys())
    else:
        return list(set(fragPool.values()))

    
def write_file(input_file,dir,maxBlocks,maxSR,asMols,minFragAtoms):

    if '.smi' in input_file:
        smiles=[line.strip().split()[0] for line in open(input_file).readlines()]
        mols=[Chem.MolFromSmiles(s) for s in smiles]
    elif '.sdf' in input_file:
        mols = Chem.SDMolSupplier(input_file)
    else:
        print ('Error: the input file should be .smi of .sdf file')
          
    if asMols == 'True':
        out_file=dir+'fragments.sdf'
        fw=Chem.SDWriter(out_file)
        for mol in mols:
            if mol is not None:
                frags=MacFrag(mol,maxBlocks=maxBlocks,maxSR=maxSR,asMols=True,minFragAtoms=minFragAtoms)
            for f in frags:
                try:
                    fw.write(f)
                except Exception:
                    continue
        
    elif asMols == 'False':
        out_file=dir+'fragments.smi'
        fw=open(out_file,'w')
        for mol in mols:
            if mol is not None:
                frags=MacFrag(mol,maxBlocks=maxBlocks,maxSR=maxSR,asMols=False,minFragAtoms=minFragAtoms)
            for f in frags:
                fw.write(f+'\n')
        fw.close()
    
  
    
def main():
    opt = parse_args()
    input_file=opt.input_file
    dir=opt.output_path
    asMols=opt.asMols
    maxBlocks=int(opt.maxBlocks)
    maxSR=int(opt.maxSR)
    minFragAtoms=int(opt.minFragAtoms)
    write_file(input_file,dir,maxBlocks,maxSR,asMols,minFragAtoms)
    
if __name__ == '__main__':    
    main()  
