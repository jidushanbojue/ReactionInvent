import pandas as pd
import pathlib
# from reaction_utils.amol_utils_rdchiral import Reaction
from reaction_utils.amol_utils_rdchiral import Reaction
from rdkit.Chem import Draw, AllChem
from multiprocessing import Pool
from rdkit import Chem

class NotCanonicalizableSmilesException(ValueError):
    pass

def canonicalize_smi(smi, remove_atom_mapping=False):
    r"""
    Canonicalize SMILES
    """
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        raise NotCanonicalizableSmilesException("Molecule not canonicalizable")
    if remove_atom_mapping:
        for atom in mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                atom.ClearProp("molAtomMapNumber")
    return Chem.MolToSmiles(mol)

def process_reaction(rxn):
    """
    Process and canonicalize reaction SMILES
    """
    reactants, reagents, products = rxn.split(">")
    try:
        precursors = [canonicalize_smi(r, True) for r in reactants.split(".")]
        if len(reagents) > 0:
            precursors += [
                canonicalize_smi(r, True) for r in reagents.split(".")
            ]
        products = [canonicalize_smi(p, True) for p in products.split(".")]
    except NotCanonicalizableSmilesException:
        return ""

    joined_precursors = ".".join(sorted(precursors))
    joined_products = ".".join(sorted(products))
    return f"{joined_precursors}>>{joined_products}"


def draw_reaction(reaction_file):
    df = pd.read_csv(reaction_file)
    for idx, line in df.iterrows():
        print(idx)

        rxn_canonical_smile = process_reaction(line['new_reaction'])
        rxn = AllChem.ReactionFromSmarts(rxn_canonical_smile)
        d2d_template = Draw.MolDraw2DCairo(800, 300)
        d2d_template.DrawReaction(rxn)
        png_template = d2d_template.GetDrawingText()
        open('./new_reaction_plot/'+str(idx)+'.png', 'wb+').write(png_template)

        # rxn_canonical_smile = process_reaction(line['rsmi'])
        # rxn_reaction_canonical = AllChem.ReactionFromSmarts(rxn_canonical_smile)
        # rxn_reaction_hasatommap = AllChem.ReactionFromSmarts(line['rsmi'])
        # # rxn_reaction = AllChem.ReactionFromSmarts(line['rsmi'])
        # rxn_template = AllChem.ReactionFromSmarts(line['retro_template'])
        #
        # d2d_reaction = Draw.MolDraw2DCairo(1000, 800)
        # d2d_reaction.DrawReaction(rxn_reaction_hasatommap)
        # png_reaction = d2d_reaction.GetDrawingText()
        # open('./reaction_plot/'+str(idx)+'_reaction_has_atommap.png', 'wb+').write(png_reaction)
        #
        # d2d_reaction = Draw.MolDraw2DCairo(1000, 800)
        # d2d_reaction.DrawReaction(rxn_reaction_canonical)
        # png_reaction = d2d_reaction.GetDrawingText()
        # open('./reaction_plot/'+str(idx)+'_reaction_canonical.png', 'wb+').write(png_reaction)
        #
        # d2d_template = Draw.MolDraw2DCairo(800, 300)
        # d2d_template.DrawReaction(rxn_template)
        # png_template = d2d_template.GetDrawingText()
        # open('./reaction_plot/'+str(idx)+'_template.png', 'wb+').write(png_template)


def test_process_info(reaction_smarts):
    reaction = Reaction(reaction_smarts)
    info = reaction.generate_reaction_template(radius=1)
    print('Done')


# def are_product_one_fragment(product, changed_atoms, changed_atom_tags):
def worker(line):
    try:
        reaction = Reaction(line['rsmi'])
        canonical_template, retro_template, changed_atoms, changed_atom_tags, reactant_fragments, product_fragments = reaction.generate_reaction_template(
            radius=1)
        if len(reactant_fragments.split('.')) == 2 and len(product_fragments.split('.')) == 1:
            return line
        else:
            return ''
    except Exception as e:
        print(e)
        return ''



def get_reaction_info(reaction_file, new_reaction_file):
    df = pd.read_csv(reaction_file)
    lines_list = []
    for idx, line in df.iterrows():
        print(idx)
        if idx < 3:
            continue
        lines_list.append(line)

    # p = Pool(160)
    # res = p.map(worker, lines_list)
    # p.close()
        # print(idx)
        reaction = Reaction(line['rsmi'])
        canonical_template, retro_template, changed_atoms, changed_atom_tags, reactant_fragments, product_fragments, leaving_group_info = reaction.generate_reaction_template(radius=1)
        if len(reactant_fragments.split('.')) == 2 and len(product_fragments.split('.')) == 1:
            lines_list.append(line)
    # res = [el for el in res if el is not '']
    # new_df = pd.DataFrame(res)
    # new_df.to_csv(new_reaction_file)

    print('Done')


def justify_NuEu(reaction_file, result_file):
    df = pd.read_csv(reaction_file)
    Eu_list = []
    Eu_frag_list = []
    Eu_changed_site_list = []
    Nu_list = []
    Nu_frag_list = []
    Nu_changed_site_list = []
    for idx, line in df.iterrows():
        print('This idx is: ', idx)
        # if idx < 194:
        #     continue
        reaction = Reaction(line['rsmi'])
        info = reaction.generate_reaction_template(radius=1)
        # frags_info = info[6]
        # print(len(frags_info))

        if info is not None and len(info[6]) == 2:
            frags_info = info[6]

            if frags_info[0]['leaving_group_atomic_num'] > frags_info[1]['leaving_group_atomic_num']:
                Eu_smiles = frags_info[0]['smiles']
                Eu_frag = frags_info[0]['frags_smiles']
                Eu_changed_tag = frags_info[0]['frag_changed_tag']

                Nu_smiles = frags_info[1]['smiles']
                Nu_frag = frags_info[1]['frags_smiles']
                Nu_changed_tag = frags_info[1]['frag_changed_tag']
                #
                # Eu_frag_list.append(Eu_frag)
                # Eu_changed_site_list.append(Eu_changed_tag)
                # Nu_frag_list.append(Nu_frag)
                # Nu_changed_site_list.append(Nu_changed_tag)

            else:
                Eu_smiles = frags_info[1]['smiles']
                Eu_frag = frags_info[1]['frags_smiles']
                Eu_changed_tag = frags_info[1]['frag_changed_tag']

                Nu_smiles = frags_info[0]['smiles']
                Nu_frag = frags_info[0]['frags_smiles']
                Nu_changed_tag = frags_info[0]['frag_changed_tag']

                # Eu_frag_list.append(Eu_frag)
                # Eu_changed_site_list.append(Eu_changed_tag)
                # Nu_frag_list.append(Nu_frag)
                # Nu_changed_site_list.append(Nu_changed_tag)

            Eu_list.append(Eu_smiles)
            Eu_frag_list.append(Eu_frag)
            Eu_changed_site_list.append(Eu_changed_tag)

            Nu_list.append(Nu_smiles)
            Nu_frag_list.append(Nu_frag)
            Nu_changed_site_list.append(Nu_changed_tag)

        else:
            Eu_list.append('')
            Eu_frag_list.append('')
            Eu_changed_site_list.append(None)
            Nu_list.append('')
            Nu_frag_list.append('')
            Nu_changed_site_list.append(None)

            # pass
            # df['Eu_frag'] = ''
            # df['Eu_changed_tag'] = None
            # df['Nu_frag'] = ''
            # df['Nu_changed_tag'] = None
    df['Eu_smiles'] = Eu_list
    df['Eu_frag'] = Eu_frag_list
    df['Eu_changed_tag'] = Eu_changed_site_list

    df['Nu_smiles'] = Nu_list
    df['Nu_frag'] = Nu_frag_list
    df['Nu_changed_tag'] = Nu_changed_site_list
    df.to_csv(result_file)

    print('Done')


    # print(df.head(5))

from rdkit.Chem import RenumberAtoms

def get_change_site_id(mol, tag):
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            if int(atom.GetProp('molAtomMapNumber')) == tag:
                changed_site_idx = atom.GetIdx()
                return changed_site_idx

def get_renumber_mol(mol, tag):
    changed_site_idx = get_change_site_id(mol, tag)
    idx_list = list(range(len(mol.GetAtoms())))
    idx_list.pop(changed_site_idx)
    idx_list = [changed_site_idx] + idx_list
    new_mol = RenumberAtoms(mol, idx_list)
    new_smiles = '*' + Chem.MolToSmiles(new_mol, canonical=False)
    return Chem.MolFromSmiles(new_smiles)


def combine_fragment(frag1_info, frag2_info):
    frag1_mol = Chem.MolFromSmiles(frag1_info[1], sanitize=False)
    frag1_tag = frag1_info[2]

    frag2_mol = Chem.MolFromSmiles(frag2_info[1], sanitize=False)
    frag2_tag = frag2_info[2]
    frag2_changed_site_id = get_change_site_id(frag2_mol, frag2_tag)

    frag1_new_mol = get_renumber_mol(frag1_mol, frag1_tag)

    rms = AllChem.ReplaceSubstructs(frag1_new_mol,
                                    Chem.MolFromSmiles('*'),
                                    frag2_mol,
                                    replacementConnectionPoint=frag2_changed_site_id)
    new_product = Chem.MolToSmiles(rms[0])

    return new_product



import itertools as it

def generate_new_reaction(src_file, new_reaction_file):
    df = pd.read_csv(src_file)
    index = df['ID']
    df_Eu = df['Eu_smiles']
    df_Eu_frag = df['Eu_frag']
    df_Eu_tag = df['Eu_changed_tag']

    df_Nu = df['Nu_smiles']
    df_Nu_frag = df['Nu_frag']
    df_Nu_tag = df['Nu_changed_tag']
    zip_Eu = list(zip(index, df_Eu_frag, df_Eu_tag, df_Eu))
    zip_Nu = list(zip(index, df_Nu_frag, df_Nu_tag, df_Nu))

    ID_list = []
    reaction_list = []
    for idx, Eu_info in enumerate(zip_Eu):
        print('This is ', idx)

        for Nu_info in zip_Nu:
            if Eu_info[0] == Nu_info[0]:
                continue

            Eu_smiles = Eu_info[3]
            Nu_smiles = Nu_info[3]
            product = combine_fragment(Eu_info, Nu_info)
            ID = Eu_info[0] + '--' + Nu_info[0]
            new_reaction = Eu_smiles + '.' + Nu_smiles + '>>' + product
            ID_list.append(ID)
            reaction_list.append(new_reaction)

    new_df = pd.DataFrame()
    new_df['ID'] = ID_list
    new_df['new_reaction'] = reaction_list
    new_df.to_csv(new_reaction_file)
    print('Done')



# print(atom.GetAtomNum())

if __name__ == '__main__':
    base_dir = pathlib.Path('./data')
    # reaction_file = pathlib.Path('./data/Cu_catalyst_727.csv')
    # new_reaction_file = pathlib.Path('./data/Cu_catalyst_727_result.csv')

    reaction_file = pathlib.Path('./data/Cu_catalyst_new_reaction.csv')
    draw_reaction(reaction_file)
    # get_reaction_info(reaction_file, new_reaction_file)
    # reaction_item = 'Br-[CH2;D2;+0:1]-[CH2;D2;+0:2]-Br.[OH;D1;+0:3]-[c:4]:[c:5]-[OH;D1;+0:6]>>[CH2;D2;+0:1]1-[CH2;D2;+0:2]-[O;H0;D2;+0:6]-[c:5]:[c:4]-[O;H0;D2;+0:3]-1'
    # test_process_info(reaction_item)
    # justify_NuEu(reaction_file, new_reaction_file)

    invent_reaction_file = pathlib.Path('./data/Cu_catalyst_665.csv')
    new_invent_reaction_file = pathlib.Path('./data/Cu_catalyst_new_reaction.csv')
    generate_new_reaction(invent_reaction_file, new_invent_reaction_file)




