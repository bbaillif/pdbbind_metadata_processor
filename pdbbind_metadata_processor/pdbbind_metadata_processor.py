import os
import pandas as pd
import urllib.request
from collections import Counter
from rdkit import Chem
from tqdm import tqdm
from Bio.PDB import PDBParser

# To be able to save conformer properties
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

class PDBIDNotInPDBBindException(Exception) :
    """Raised if the input PDB ID is not in PDBBind
    
    :param pdb_id: Input PDB ID
    :type pdb_id: str
    """
    
    def __init__(self, pdb_id):
        message = f'{pdb_id} is not in PDBBind'
        super().__init__(message)


class PDBBindMetadataProcessor() :
    """Object to handle PDBBind metadata
    :param root: directory where PDBBind data are located
    :type root: str
    """
    
    def __init__(self, 
                 root: str='/home/benoit/PDBBind/',
                 corrected: bool=False) :
        self.root = root
        self.corrected = corrected
        
        self.general_dir_path = os.path.join(self.root, 
                                            'PDBbind_v2020_other_PL',
                                            'v2020-other-PL/')
        self.refined_dir_path = os.path.join(self.root, 
                                            'PDBbind_v2020_refined',
                                            'refined-set/')
        
        if self.corrected :
            self.pdb_dir = os.path.join(self.root,
                                        'new_pdbs')
            if not os.path.exists(self.pdb_dir) :
                os.mkdir(self.pdb_dir)
            self.ligand_dir = os.path.join(self.root,
                                        'new_ligands')
            if not os.path.exists(self.ligand_dir) :
                os.mkdir(self.ligand_dir)
            self.available_files = os.listdir(self.ligand_dir)
            self.available_structures = [file.split('.mol2')[0]
                                         for file in self.available_files]
        
        else :
            self.general_available_structures = os.listdir(self.general_dir_path)
            self.refined_available_structures = os.listdir(self.refined_dir_path)
            self.available_structures = set(self.general_available_structures
                                            + self.refined_available_structures)
            self.available_structures.remove('index')
        #self.chembl_targets_df = pd.read_csv('chembl_targets.csv', sep=';')
        
    def get_master_dataframe(self, remove_peptide_ligands=True) :
        """Creates a single dataframe compiling all available metadata
        in PDBBind (PL_data and PL_name files)
        :param remove_peptide_ligands: if True, remove the ligand whose name
            includes "mer"
        :type remove_peptide_ligands: bool
        :return: A single pandas DataFrame compiling metadata
        :rtype: pandas.DataFrame
        """
        
        pl_data = self.get_pl_data()
        pl_name = self.get_pl_name()
        
        pl_data['activity_list'] = pl_data['Kd/Ki'].apply(self.parse_activity_string)
        pl_data['sep'] = pl_data['activity_list'].apply(lambda x : x[1])
        pl_data['value'] = pl_data['activity_list'].apply(self.get_nanomolar_activity)
        pl_data['units'] = 'nM'
        pl_all = pl_data.merge(pl_name, on='PDB code')
        active_threshold = 1000 # 1 uM = 1000 nM
        pl_all['active'] = ((pl_all['value'] < active_threshold) 
                            & ~(pl_all['sep'].isin(['>', '~'])))
        
        if remove_peptide_ligands :
            pl_all = pl_all[~(pl_all['ligand name'].str.contains('mer'))]
            
        self.pl_all = pl_all
        return pl_all
        
    def get_pl_data(self) :
        """Compile the PL (protein-ligand) data from PL_data file
        
        :return: A single pandas DataFrame formatting the PL_data file
        :rtype: pandas.DataFrame
        """
        widths = [6,6,7,6,17,9,200]
        cols = ['PDB code',
                'resolution',
                'release year',
                '-logKd/Ki',
                'Kd/Ki',
                'reference',
                'ligand name']
        file_path = os.path.join(self.general_dir_path, 
                                 'index', 
                                 'INDEX_general_PL_data.2020')
        pl_data = pd.read_fwf(file_path, widths=widths, skiprows=6, header=None)
        pl_data.columns=cols
        return pl_data
    
    
    def get_pl_name(self) :
        """Compile the PL (protein-ligand) data from PL_name file
        
        :return: A single pandas DataFrame formatting the PL_name file
        :rtype: pandas.DataFrame
        """
        widths = [6,6,8,200]
        cols = ['PDB code',
                'release year',
                'Uniprot ID',
                'protein name']
        file_path = os.path.join(self.general_dir_path, 
                                 'index', 
                                 'INDEX_general_PL_name.2020')
        pl_name = pd.read_fwf(file_path, widths=widths, skiprows=6, header=None)
        pl_name.columns=cols
        return pl_name
    
    
    def get_pl_general(self) :
        """Compile the PL (protein-ligand) data from PL_general file
        
        :return: A single pandas DataFrame formatting the PL_general file
        :rtype: pandas.DataFrame
        """
        widths = [6,6,6,17,9,200]
        cols = ['PDB code', 
                'resolution',
                'release year', 
                'binding data', 
                'reference',
                'ligand name']
        file_path = os.path.join(self.general_dir_path, 
                                 'index', 
                                 'INDEX_general_PL.2020')
        pl_general = pd.read_fwf(file_path, widths=widths,skiprows=6,header=None)
        pl_general.columns=cols
        return pl_general
    
    
    def find_sep_in_activity_string(self, 
                                    string, 
                                    possible_seps=['=', '<', '>']) :
        """In PDBBind, the activity is represented by a string in Kd/Ki column. 
        This function finds the separator which is the sign in the 
        activity equality.
        
        :param string: Input string from the activity field in master table.
        :type string: str
        :param possible_seps: Possible separators to find in the string. 
            Defaults are "=", "<" or ">".
        :type possible_seps: list[str]
        :return: Found separator. Default is "~" if it does not find any of the
            possible separators
        :rtype: str
        """
        found_sep = '~' # default value
        for sep in possible_seps :
            if sep in string :
                found_sep = sep
                break
        return found_sep

    
    def parse_activity_string(self, string) :
        """Parse the activity (Kd/Ki) string
        Example : "Ki=400mM //" should return ["mM", "=", "400"]
        
        :param string: Activity string to parse
        :type string: str
        :return: A list [unit, sep, value]
        :rtype: list[str, str, str]
        """
        # TODO: Add the type of activity in the return(Kd or Ki or ?)
        
        sep = self.find_sep_in_activity_string(string)
        splitted_string = string.split(sep)
        value_unit = splitted_string[1]

        # maybe a better way to do this with for loop
        parsed = False
        i = 0
        value = ''
        units = ''
        while not parsed :
            char = value_unit[i]
            if char in '0123456789.' :
                value = value + char
            elif char == ' ' :
                parsed = True
            else :
                units = units + char
            i = i + 1
        return [units, sep, value]

    
    def get_nanomolar_activity(self, l) :
        """Get all activity in nanomolar
        
        :param l: list [unit, sep, value] representing activity (i.e. output 
        from parse_activity_string function)
        :type l: list[str, str, str]
        :return: Activity value in nanomolar (nM)
        :rtype: float
        """
        
        unit, sep, value = l
        value = float(value)
        
        # can be done with match case since 3.10
        if unit == 'uM' :
            value = value * 1000
        elif unit == 'mM' :
            value = value * 1000000
        elif unit == 'pM' :
            value = value / 1000
        elif unit == 'fM' :
            value = value / 1000000

        return value

    
    def get_training_test_sets(self, 
                               mode='ligand',
                               n_classes=50, 
                               train_ratio=0.6) :
        """Performs a train test split per ligand or protein based on the
        complementary binder. For a given ligand/protein name, the list of 
        protein/ligand it binds are stored, reverse sorted by number of 
        occurences, and filling first the training set and then the test set.
        
        :param n_classes: Number of ligand to include in the dataset 
            (default 50)
        :type n_classes: int
        :param train_ratio: Minimum ratio of samples in the training set
        :type train_ratio: float
        :return: Two dicts train_set and test_set, storing for each ligand key
            a list of pdb_ids
        :rtype: tuple(dict[list], dict[list])
        """
        assert mode in ['ligand', 'protein']
        if mode == 'ligand' :
            class_column_name = 'ligand name'
            binder_column_name = 'Uniprot ID'
        elif mode == 'protein' :
            class_column_name = 'Uniprot ID'
            binder_column_name = 'ligand name'
        
        train_set = {}
        test_set = {}
        pl_all = self.pl_all
        if mode == 'protein' :
            pl_all = pl_all[pl_all['Uniprot ID'] != '------']
        class_counts = pl_all[class_column_name].value_counts()
        topN_class_counts = class_counts[:n_classes]
        for class_name in topN_class_counts.index :
            pl_class = pl_all[pl_all[class_column_name] == class_name]
            pl_class = pl_class[pl_class['PDB code'].isin(self.available_structures)]

            # Make sure we have enough data for given class
            if len(pl_class) > 10 : 
                train_pdbs = []
                test_pdbs = []
                counter = Counter()
                counter.update(pl_class[binder_column_name].values)
                if len(counter) > 1 :
                    for binder_name, count in counter.most_common() :
                        pdb_ids = pl_class[pl_class[binder_column_name] == binder_name]['PDB code'].values
                        if len(train_pdbs) < len(pl_class) * train_ratio :
                            train_pdbs.extend(pdb_ids)
                        else :
                            test_pdbs.extend(pdb_ids)
                    train_set[class_name] = train_pdbs
                    test_set[class_name] = test_pdbs
                    
        return (train_set, test_set)
    
    
    def get_pdb_id_pathes(self, 
                          pdb_id: str, 
                          ligand_format: str='sdf') :
        """Give the path to the protein pdb and ligand sdf file(s) for a
        given pdb_id if present in PDBbind
        
        :param pdb_id: Input PDB ID
        :type pdb_id: str
        :param ligand_format: Format of the ligand to return (sdf or mol2)
        :type ligand_format: str
        :return: Tuple with the protein path and the ligand path(es)
        :rtype: tuple(str, list[str])
        """
        
        assert ligand_format in ['sdf', 'mol2'], 'Ligand format is sdf or mol2'
        
        if self.corrected :
            
            if pdb_id in self.available_structures :
                protein_path = os.path.join(self.pdb_dir, 
                                            f'{pdb_id}.pdb')
                ligand_path = os.path.join(self.ligand_dir, 
                                           f'{pdb_id}.mol2')
            else :
                raise PDBIDNotInPDBBindException(pdb_id)
            
        else :
            if pdb_id in self.general_available_structures :
                correct_dir_path = self.general_dir_path
            elif pdb_id in self.refined_available_structures :
                correct_dir_path = self.refined_dir_path
            else :
                raise PDBIDNotInPDBBindException(pdb_id)
            
            protein_path = os.path.join(correct_dir_path, 
                                        pdb_id, 
                                        f'{pdb_id}_protein.pdb')
            ligand_path = os.path.join(correct_dir_path, 
                                        pdb_id, 
                                        f'{pdb_id}_ligand.{ligand_format}')
            
        ligand_pathes = [ligand_path]
        return protein_path, ligand_pathes
    
    
    def get_ligand_name(self, 
                        pdb_id: str) :
        master_table = self.get_master_dataframe()
        pdb_line = master_table[master_table['PDB code'] == pdb_id]
        ligand_name = pdb_line['ligand name'].values[0]
        return ligand_name
    
    
    def get_protein_name(self,
                         uniprot_id: str) :
        master_table = self.get_master_dataframe()
        pdb_lines = master_table[master_table['Uniprot ID'] == uniprot_id]
        protein_names = pdb_lines['protein name'].values
        counter = Counter(protein_names)
        return counter.most_common()[0][0]
    
    def get_chains(self,
                   pdb_id: str) :
        protein_path, ligand_pathes = self.get_pdb_id_pathes(pdb_id)
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure('struct', protein_path)
        chains = []
        for model in structure :
            for chain in model :
                chains.append(chain.id)
        return set(chains)
    
    def get_molecules(self, 
                      subset='all',
                      filter_mers=False,
                      ligand_format='mol2') :
        
        assert subset in ['all', 'general', 'refined']
        if subset == 'all' :
            pdb_ids = self.available_structures
        elif subset == 'general' :
            pdb_ids = self.general_available_structures
        elif subset == 'refined' :
            pdb_ids = self.refined_available_structures
        
        if filter_mers :
            table = self.get_master_dataframe()
            pdb_ids = [pdb_id 
                       for pdb_id in pdb_ids 
                       if pdb_id in table['PDB code'].values]
        
        mols = []
        for pdb_id in tqdm(pdb_ids) :
            protein_path, ligand_pathes = self.get_pdb_id_pathes(pdb_id=pdb_id,
                                                                 ligand_format=ligand_format)
            ligand_path = ligand_pathes[0]
            try :
                with open(ligand_path, 'r') as f :
                    mol2block = f.readlines()
                #import pdb;pdb.set_trace()
                mol2block = [line.replace('CL', 'Cl') for line in mol2block]
                mol2block = [line.replace('BR', 'Br') for line in mol2block]
                mol2block = [line.replace('AS', 'As') for line in mol2block]
                mol2block = ''.join(mol2block)
                mol = Chem.rdmolfiles.MolFromMol2Block(mol2block)
                if mol is not None :
                    rdmol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
                    if rdmol is not None : #rdkit parsable
                        mol.GetConformer().SetProp('PDB_ID', pdb_id)
                        mol.GetConformer().SetProp('pdbbind_id', pdb_id)
                        mols.append(mol)
                    else :
                        print(f'{pdb_id} Not RDKit parsable')
            except Exception as e :
                print('Impossible to read mol2 file for ' + pdb_id)
                print(str(e))
        return mols
    
    def download_all_corrected_data(self) :
        table = self.get_master_dataframe()
        for i, row in table.iterrows() :
            pdb_id = row['PDB code']
            ligand_name = row['ligand name']
            ligand_name = ligand_name[1:-1] # remove parenthesis
            self.download_corrected_data(pdb_id,
                                         ligand_name)
    
    
    def download_corrected_data(self, 
                                pdb_id,
                                ligand_name) :
        self.download_pdb_file(pdb_id)
        pdb_path = os.path.join(self.pdb_dir,
                                f'{pdb_id}.pdb')
        with open(pdb_path, 'r') as f :
            lines = f.readlines()
        for line in lines :
            if line.startswith('HETATM') :
                ligand_entry = line[17:20].strip()
                chain_entry = line[20:22].strip()
                res_entry = line[22:26].strip()
                if ligand_entry == ligand_name :
                    self.download_ligand_file(pdb_id, 
                                            chain=chain_entry, 
                                            res_id=res_entry)
                    break

        
    def download_pdb_file(self, 
                          pdb_id) :
        pdb_path = os.path.join(self.pdb_dir,
                                f'{pdb_id}.pdb')
        if not os.path.exists(pdb_path) :
            urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdb_id}.pdb', 
                                       pdb_path)
            
            
    def download_ligand_file(self,
                             pdb_id,
                             chain,
                             res_id) :
        ligand_path = os.path.join(self.ligand_dir,
                                   f'{pdb_id}.mol2')
        if not os.path.exists(ligand_path) :
            url = f'https://models.rcsb.org/v1/{pdb_id}/ligand?auth_asym_id={chain}&auth_seq_id={res_id}&encoding=mol2'
            #import pdb;pdb.set_trace()
            urllib.request.urlretrieve(url, 
                                       ligand_path)
            
            
    def generate_complex(self,
                         pdb_id,
                         force_write=True) :
        protein_path, ligand_pathes = self.get_pdb_id_pathes(pdb_id=pdb_id)
        ligand_path = ligand_pathes[0]
        complex_path = self.get_complex_path(pdb_id=pdb_id)
        if not os.path.exists(complex_path) or force_write :
            from pymol import cmd
            cmd.set('retain_order', 1)
            cmd.load(protein_path)
            cmd.load(ligand_path)
            cmd.save(complex_path)
            cmd.delete('all')
        
        
    def generate_complexes(self) :
        for pdb_id in tqdm(self.available_structures) :
            self.generate_complex(pdb_id=pdb_id)
            
            
    def get_complex_path(self, pdb_id) :
        if pdb_id in self.general_available_structures :
            correct_dir_path = self.general_dir_path
        elif pdb_id in self.refined_available_structures :
            correct_dir_path = self.refined_dir_path
        else :
            raise PDBIDNotInPDBBindException(pdb_id)
        
        complex_path = os.path.join(correct_dir_path, 
                                    pdb_id, 
                                    f'{pdb_id}_complex.pdb')
        
        return complex_path