# Function to lookup the targets of a given compound
# in the ChEMBL database using their REST API. Returns
# a list of gene symbols. 
def query_chembl(inchi):
    import requests
    human_target_genes = []
    inchi_url = "https://www.ebi.ac.uk/chembl/api/data/similarity/{}/100?format=json".format(inchi)
    inchi_request_json = requests.get(inchi_url).json()
    if 'molecules' in inchi_request_json.keys(): # if the compound is in the CHEMBL database
        chembl_ids = [ m['molecule_chembl_id'] for m in inchi_request_json['molecules'] ]

        # see if the CHEMBL_IDs correspond to any known mechanisms of action (targets)
        mech_url_temp = "https://www.ebi.ac.uk/chembl/api/data/mechanism?molecule_chembl_id={}&format=json"
        for cid in chembl_ids:
            mech_url = mech_url_temp.format(cid)
            mech_request_json = requests.get(mech_url).json()
            if mech_request_json['mechanisms']:

                # if there are known mechanisms (targets), find their gene symbols
                for mech in mech_request_json['mechanisms']:
                    target_cid = mech['target_chembl_id']
                    target_url = "https://www.ebi.ac.uk/chembl/api/data/target?target_chembl_id={}&format=json".format(target_cid)
                    target_request_json = requests.get(target_url).json()

                    # filter out non-human targets
                    human_targets = [ t for t in target_request_json['targets'] if t['organism'] == 'Homo sapiens' ]
                    # targets can be protein complexes, break them down into individual components (genes)
                    human_target_genes += [ htcs['component_synonym'] 
                                               for ht in human_targets 
                                               for htc in  ht['target_components'] 
                                               for htcs in htc['target_component_synonyms']
                                               if htcs['syn_type'] == 'GENE_SYMBOL'
                                          ]
    return human_target_genes                


# Function to parse DrugBank database XML file
# and convert it to a dictionary for drug target lookup.
# Returns a dictionary whose keys are the drugs' InChiKeys 
# and whose values are lists of the drugs' targets' Hugo Gene Symbols.
def clean_drugbank(db_path):
    import xmltodict
    drugbank_db = xmltodict.parse(open(db_path).read())
    drug0 = drugbank_db['drugbank']['drug'][0] # an example drug to use as a type reference
    drugbank_target_dict = {}

    # all drugs are OrderedDicts 
    for drug in drugbank_db['drugbank']['drug']:
        
        # limit analysis to small molecules (n=9292), ignore type 'biotech' (n=1213)
        if drug['@type'] == 'small molecule':
            calc_props = drug['calculated-properties']
            
            # find the InchiKey for the drug
            # n=576 cpds have 'NoneType' props
            # the rest are OrderedDicts w/ exactly one key: 'property'
            if calc_props: 
                props = calc_props['property'] 
                # props is usually a list of OrderedDicts, each with keys 'kind', 'value', and 'source'
                # DB11635 is an exception, for which props is not a list, but a single OrderedDict with no Inchi
                if type(props) == type([]):
                    inchi_list = [ prop['value'] for prop in props if prop['kind'] == 'InChIKey' ]
                    # all cpds have either 1 or 0 inchi_keys
                    if len(inchi_list) > 0:
                        inchi_key = inchi_list[0]
                        gene_targets = []
                        
                        # Now we have to lookup the drug targets. 
                        # All drugs have a 'target' field, but n=2259 have are NoneType
                        targets = drug['targets']
                        if targets: # (n=6455) targets is an OrderedDict with one key: 'target'
                            target = targets['target'] # can either be an OrderedDict or a list of OrderedDicts
                            
                            # target is an OrderedDict (n=4152)
                            if type(target) == type(drug0): 
                                # some targets are not proteins (eg: DNA, RNA, ribosome)
                                # we'll only pay attention to protein targets 
                                if 'polypeptide' in target.keys():
                                    polypeptide = target['polypeptide'] # polypeptide can be either a list or an OrderedDict
                                    
                                    # polypeptide is an OrderedDict
                                    if type(polypeptide) == type(drug0):
                                        gene_name =  polypeptide['gene-name']
                                        if gene_name and gene_name.isupper(): # restrict to official gene symbols (i.e. not 'pol' or 'HIV protease')
                                            gene_targets.append(gene_name)
                                        
                                    # polypeptide is a list of OrderedDicts (n=14)
                                    elif type(polypeptide) == type([]):
                                        gene_names = [ poly['gene-name'] for poly in polypeptide ]
                                        gene_targets += gene_names                                    
                                    
                            # target is a list of OrderedDicts (n=2303)
                            elif type(target) == type([]):
                                for targ in target:
                                    if 'polypeptide' in targ.keys():
                                        polypeptide = targ['polypeptide'] # polypeptide can be either a list or an OrderedDict
                                        
                                        # polypeptide is an OrderedDict
                                        if type(polypeptide) == type(drug0):
                                            gene_name =  polypeptide['gene-name']
                                            if gene_name and gene_name.isupper(): # restrict to official gene symbols (i.e. not 'pol' or 'HIV protease')
                                                gene_targets.append(gene_name)

                                        # polypeptide is a list of OrderedDicts (n=135)
                                        elif type(polypeptide) == type([]):
                                            gene_names = [ poly['gene-name'] for poly in polypeptide ]
                                            gene_targets += gene_names          
                        
                        if len(gene_targets) > 0:
                            drugbank_target_dict[inchi_key] = gene_targets

    return drugbank_target_dict


# Widget based progress bar for Jupyter (IPython Notebook)
# https://github.com/alexanderkuk/log-progress
def log_progress(sequence, every=None, size=None, name='Items'):
    from ipywidgets import IntProgress, HTML, VBox
    from IPython.display import display

    is_iterator = False
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = int(size / 200)     # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    display(box)

    index = 0
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{name}: {index} / ?'.format(
                        name=name,
                        index=index
                    )
                else:
                    progress.value = index
                    label.value = u'{name}: {index} / {size}'.format(
                        name=name,
                        index=index,
                        size=size
                    )
            yield record
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'success'
        progress.value = index
        label.value = "{name}: {index}".format(
            name=name,
            index=str(index or '?')
        )


