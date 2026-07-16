#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas
from PIL import ImageDraw
from PIL import ImageFont
import os
import sys

if len(sys.argv) == 2:
    smirks_fn = sys.argv[1]
else:
    smirks_fn = "smirks.txt"

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# old code path
from tautomerizer import Tautomerizer as Tautomerizer2021
tautomerizer = Tautomerizer2021(smirks_fn)
old_tautomerizer = tautomerizer.get_tautomers

# from modern molscrub
from molscrub import Tautomerizer

#font = ImageFont.truetype("FreeMono.ttf", 16)
#font_small = ImageFont.truetype("FreeMono.ttf", 16)

df = pandas.read_csv("Tautomer_database_release_3a.csv")

imgsize = 300


largest_Fragment = rdMolStandardize.LargestFragmentChooser()


def evaluate_row(tautomers_list_of_lists, mols, is_wanted):
    target_smiles_list = []
    for i, mol in enumerate(mols):
        if not is_wanted[i]:
            continue
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)), isomericSmiles=False)
        target_smiles_list.append(smiles)
    n_hits = 0
    for tautomers in tautomers_list_of_lists:
        for pred_smiles in tautomers:
            n_hits += int(pred_smiles in target_smiles_list)
    return n_hits

def tautomerize_row(mols, tautomer_generator_func):
    tautomers_list_of_lists = []
    for mol in mols:
        tautomers = tautomer_generator_func(mol)
        smiles_set = set()
        for tautomer in tautomers: 
            tautomer = Chem.RemoveHs(tautomer)
            pred_smiles = Chem.MolToSmiles(tautomer, isomericSmiles=False)
            smiles_set.add(pred_smiles)
        tautomers_list_of_lists.append(smiles_set)
    return tautomers_list_of_lists


def evalondf(df, wanted_transforms, tautomerizer, category_thresholds):

    stats = {
        "generated": 0,
        "wanted": {t: 0 for t in category_thresholds},
        "hits": {t: 0 for t in category_thresholds},
        "rows_skipped": 0,
        "tautomers_skipped": 0,
        "tautomers_skipped_in_kept_rows": 0,
    }

    entries_with_none = 0
    for index, row  in df.iterrows():
        #if index!=689:continue
        tautomers = [Chem.MolFromSmiles(row["SMILES_%d" % (i+1)]) for i in range(row["Size"])]
        if None in tautomers:
            print("Mol is None, index=%d, skipping" % index)
            entries_with_none += 1
            continue
    
        all_tautomers = [largest_Fragment.choose(mol) for mol in tautomers]
        tautomers = [all_tautomers[0]]
        kept_tautomer_indices = [0]
        skipped_tautomers = 0
        for i in range(row["Size"]):
            #k = row["Quantitative_ratio_%d" % (i+1)]
            #q = row["Qualitative_prevalence_%d" % (i+1)]
            if i > 0:
                transform = row["Transf_1_%d" % (i+1)]
                if wanted_transforms == "all" or transform in wanted_transforms:
                    tautomers.append(all_tautomers[i])
                    kept_tautomer_indices.append(i)
                else:
                    skipped_tautomers += 1

        is_wanted = {t: [] for t in category_thresholds}
        for i in kept_tautomer_indices:
            prevalence = row["Prevalence_Category_%d" % (i+1)]
            for threshold in category_thresholds:
                is_wanted[threshold].append(prevalence >= threshold)


        stats["tautomers_skipped"] += skipped_tautomers
    
        if len(tautomers) == 0:
            stats["rows_skipped"] += 1
            continue
        stats["tautomers_skipped_in_kept_rows"] += skipped_tautomers
    
        try:
            tautomers_list_of_lists = tautomerize_row(tautomers, tautomerizer) 
            n_total = sum([len(t) for t in tautomers_list_of_lists])
            n_hits = {}
            n_wanted = {}
            no_bueno = False
            for threshold in category_thresholds:
                n_hits[threshold] = evaluate_row(tautomers_list_of_lists, tautomers, is_wanted[threshold])
                n_wanted[threshold] = sum(is_wanted[threshold]) * len(tautomers)
                if n_wanted[threshold] < n_hits[threshold]:
                    no_bueno = True
                    for tautomer in tautomers:
                        print(Chem.MolToSmiles(tautomer))
                    print(f"{is_wanted=}")
                    print(f"{n_wanted=} {n_hits=} {n_total=} {len(tautomers)=}")
                    for ts in tautomers_list_of_lists:
                        print(len(ts))
                        all_smi = set()
                        for t in ts:
                            pred_smiles = Chem.MolToSmiles(Chem.RemoveHs(t), isomericSmiles=False)
                            all_smi.add(pred_smiles)
                            print(f"{pred_smiles=}")
                        print(f"{len(all_smi)=}")
        
        except:
            n_hits = {t: 0 for t in category_thresholds}
            n_wanted = {t: 0 for t in category_thresholds}
            n_total = 0
            print(f"oops failed {index=}")
        if no_bueno:
            print("Exiting because nr hits exceeded nr of expected tautomers - something went terribly wrong.")
            sys.exit()

        stats["generated"] += n_total
        for threshold in category_thresholds:
            stats["hits"][threshold] += n_hits[threshold]
            stats["wanted"][threshold] += n_wanted[threshold]

    return stats
    
def get_rdkit_tautomers(mol):
    """ based on: https://gist.github.com/iwatobipen/ca1999b6e4637daf88f315b412220737
    """
    tenum = rdMolStandardize.TautomerEnumerator()
    tenum.Canonicalize(mol)
    res = tenum.Enumerate(mol)
    return list(res)

def print_stats(stats, do_rdkit=False):
    string = ""
    string += "success: %d out of %d\n" % (stats['success_count'], stats['transform_count'])
    string += "failure: %d out of %d\n" % (stats['failure_count'], stats['transform_count'])
    string += "good hits:   %4d\n" % (stats['good_hit'])
    string += " bad hits:   %4d\n" % (stats['bad_hit'])
    string += "good misses: %4d\n" % (stats['good_miss'])
    string += " bad misses: %4d\n" % (stats['bad_miss'])
    string += "total generated tautomers: %d\n" % (stats['generated_count'])
    string += "DB entries with 'None' molecules (were excluded): %d\n" % (stats['entries_with_none'])
    if do_rdkit:
        string += "(rdkit) success: %d out of %d\n" % (stats['rdkit_success_count'], stats['transform_count'])
        string += "(rdkit) total generated tautomers: %d\n" % (stats['rdkit_generated_count'])
    return string

target_transforms_list = [
    "all",
    ["PT_06_00", "PT_07_00", "PT_09_00"],
    ["PT_09_00"],
    ["PT_06_00"],
    ["PT_07_00"],
]
results = {}
write_figs_misses_excesses = True

molscrub_tautomerizer = Tautomerizer.from_default_data_files()

rules_cfbf0b1 = Tautomerizer.from_reactions_filename("tautomers-cfbf0b1.txt")
rules_f714c5a = Tautomerizer.from_reactions_filename("tautomers-f714c5a.txt")
rules_fd70a70 = Tautomerizer.from_reactions_filename("tautomers-fd70a70.txt")

funcs = {
    "installed": molscrub_tautomerizer,
    "cfbf0b1": rules_cfbf0b1,
    "f714c5a": rules_f714c5a,
    "fd70a70": rules_fd70a70,
    "2021": old_tautomerizer,
    "rdkit": get_rdkit_tautomers,
}
category_thresholds = [2, 3, 4]


for target_transforms in target_transforms_list:
    for func_name, func in funcs.items():
        stats = evalondf(df, target_transforms, func, category_thresholds)
        for threshold in category_thresholds:
            hits = stats["hits"][threshold]
            wanted = stats["wanted"][threshold] 
            recall = hits / wanted
            precision = hits / stats["generated"]
            print(f"{target_transforms=} {func_name=} {threshold=} {hits=} {wanted=} {stats['generated']=} {precision=:.3f} {recall=:.3f}")
