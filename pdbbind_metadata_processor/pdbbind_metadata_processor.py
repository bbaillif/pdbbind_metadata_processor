import os
import pandas as pd
from collections import Counter

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
    
    def __init__(self, root='/home/benoit/PDBBind/') :
        self.root = root
        self.general_dir_path = os.path.join(self.root, 
                                        'PDBbind_v2020_other_PL',
                                        'v2020-other-PL/')
        self.refined_dir_path = os.path.join(self.root, 
                                        'PDBbind_v2020_refined',
                                        'refined-set/')#
        self.general_available_structures = os.listdir(self.general_dir_path)
        self.refined_available_structures = os.listdir(self.refined_dir_path)
        self.all_available_structures = (self.general_available_structures
                                         + self.refined_available_structures)
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

    
    def get_training_test_sets(self, n_ligands=50, train_ratio=0.7) :
        """Performs a train test split per ligand based on protein which it 
        binds. For a given ligand name, the list of protein it binds are stored,
        reverse sorted by number of occurences, and filling first the 
        training set and then the test set.
        
        :param n_ligands: Number of ligand to include in the dataset 
            (default 50)
        :type n_ligands: int
        :param train_ratio: Minimum ratio of samples in the training set
        :type train_ratio: float
        :return: Two dicts train_set and test_set, storing for each ligand key
            a list of pdb_ids
        :rtype: tuple(dict[list], dict[list])
        """
        train_set = {}
        test_set = {}
        ligand_counts = self.pl_all['ligand name'].value_counts()
        topN_ligand_counts = ligand_counts[:n_ligands]
        for ligand_name in topN_ligand_counts.index :
            pl_lig = self.pl_all[self.pl_all['ligand name'] == ligand_name]
            pl_lig = pl_lig[pl_lig['PDB code'].isin(self.all_available_structures)]

            # Make sure we have enough data for given ligand
            if len(pl_lig) > 10 : 
                train_pdbs = []
                test_pdbs = []
                counter = Counter()
                counter.update(pl_lig['protein name'].values)
                if len(counter) > 1 :
                    for prot_name, count in counter.most_common() :
                        pdb_ids = pl_lig[pl_lig['protein name'] == prot_name]['PDB code'].values
                        if len(train_pdbs) < len(pl_lig) * train_ratio :
                            train_pdbs.extend(pdb_ids)
                        else :
                            test_pdbs.extend(pdb_ids)
                    train_set[ligand_name] = train_pdbs
                    test_set[ligand_name] = test_pdbs
                    
        return (train_set, test_set)
    
    
    def get_pdb_id_pathes(self, pdb_id) :
        """Give the path to the protein pdb and ligand sdf files for a
        given pdb_id if present in PDBbind
        
        :param pdb_id: Input PDB ID
        :type pdb_id: str
        :return: Tuple with the protein path and the ligand path
        :rtype: tuple(str, str)
        """
        
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
                                    f'{pdb_id}_ligand.sdf')
            
        return protein_path, ligand_path