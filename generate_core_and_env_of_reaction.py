import pandas as pd
import pathlib
# from reaction_utils.amol_utils_rdchiral import Reaction
from reaction_utils.amol_utils_rdchiral import Reaction
from rdkit.Chem import Draw, AllChem
from multiprocessing import Pool

def draw_reaction(reaction_file):
    df = pd.read_csv(reaction_file)
    for idx, line in df.iterrows():
        print(idx)
        rxn_reaction = AllChem.ReactionFromSmarts(line['rsmi'])
        rxn_template = AllChem.ReactionFromSmarts(line['retro_template'])
        d2d_reaction = Draw.MolDraw2DCairo(800, 300)
        d2d_reaction.DrawReaction(rxn_reaction)
        png_reaction = d2d_reaction.GetDrawingText()
        open('./reaction_plot/'+str(idx)+'_reaction.png', 'wb+').write(png_reaction)

        d2d_template = Draw.MolDraw2DCairo(800, 300)
        d2d_template.DrawReaction(rxn_template)
        png_template = d2d_template.GetDrawingText()
        open('./reaction_plot/'+str(idx)+'_template.png', 'wb+').write(png_template)


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
        lines_list.append(line)

    p = Pool(160)
    res = p.map(worker, lines_list)
    p.close()
        # print(idx)
        # reaction = Reaction(line['rsmi'])
        # canonical_template, retro_template, changed_atoms, changed_atom_tags, reactant_fragments, product_fragments = reaction.generate_reaction_template(radius=1)
        # if len(reactant_fragments.split('.')) == 2 and len(product_fragments.split('.')) == 1:
        #     lines_list.append(line)
    res = [el for el in res if el is not '']
    new_df = pd.DataFrame(res)
    new_df.to_csv(new_reaction_file)

    print('Done')

if __name__ == '__main__':
    base_dir = pathlib.Path('./data')
    reaction_file = pathlib.Path('./data/746w_reaction.csv')
    new_reaction_file = pathlib.Path('./data/746w_reaction_result.csv')
    # draw_reaction(reaction_file)
    get_reaction_info(reaction_file, new_reaction_file)
    # reaction_item = 'Br-[CH2;D2;+0:1]-[CH2;D2;+0:2]-Br.[OH;D1;+0:3]-[c:4]:[c:5]-[OH;D1;+0:6]>>[CH2;D2;+0:1]1-[CH2;D2;+0:2]-[O;H0;D2;+0:6]-[c:5]:[c:4]-[O;H0;D2;+0:3]-1'
    # test_process_info(reaction_item)


