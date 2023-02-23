PASER_PLOT_HELP_MSG = '''   
   
    **Peptide Uniqueness**
     
    * is defined by its **(sequence, charge)** pair. For example: (PEPTIDE, 2) != (PEPTIDE, 3) != (PEPT(15.2)IDE, 2) != (PEMMIDE, 2)
    
    **Data Frame Columns**
    
    * **name**: This element represents a label associated with the experiment.
    
    * **order**: This element represents the order of the experiment.
    
    * **unique_peptides**: This element represents the number of unique peptides in the experiment.
    
    * **duplicate_peptides**: This element represents the number of duplicate peptides in the experiment, i.e., the number of peptides that appear more than once in the experiment.
    
    * **total_peptides**: This element represents the total number of peptides in the experiment.
        
    * **new_unique_peptides**: This element represents the number of new unique peptides in the experiment that are not present in the previous experiments.
    
    * **new_peptides**: This element represents the number of new peptides in the experiment that are not present in the previous experiments.
    
    * **seen_unique_peptides**: This element represents the number of unique peptides in the experiment that are also present in the previous experiments.
    
    * **seen_peptides**: This element represents the number of peptides in the experiment that are also present in the previous experiments.
    
    * **unique_proteins**: This element represents the number of unique proteins in the experiment.
    
    * **duplicate_proteins**: This element represents the number of duplicate proteins in the data, i.e., the number of proteins that appear more than once in the experiment.
    
    * **total_proteins**: This element represents the total number of proteins in the experiment.
        
    * **new_proteins**: This element represents the number of new proteins in the experiment that are not present in the previous experiments.
    
    * **new_unique_proteins**: This element represents the number of new unique proteins in the experiment that are not present in the previous experiments.
    
    * **seen_proteins**: This element represents the number of proteins in the experiment that are also present in the previous experiments.
    
    * **seen_unique_proteins**: This element represents the number of unique proteins in the experiment that are also present in the previous experiments.
    '''

PASER_VENN_HELP_MSG = '''
    
    Upload 2-3 DTASelect-filter.txt files.

    Protein Counts - number of unique protein locus's

    Peptide Counts - number of unique peptide sequences
    '''