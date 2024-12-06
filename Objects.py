

class Object():
    def __init__(self,
                 species='human',
                 gene=None,
                 ):
        
        import os
        self.fpath = os.path.join(self.dir['base'],self.dir['data'],'RefSeq',species,gene.upper(),'data_report.jsonl')
        
        import json
        with open(self.fpath, 'r') as f:
            data = json.load(f)
        self.chromosome = data['chromosomes']
        self.position = [int(data['annotations'][0]['genomicLocations'][0]['genomicRange']['begin']),\
                         int(data['annotations'][0]['genomicLocations'][0]['genomicRange']['end'])]


        self.dna = Dna()
        self.rna = Rna()
        self.protein = Protein()

class Dna(Object):
    def __init__(self):
        pass

class Rna(Object):
    def __init__(self):
        pass

class Protein(Object):
    def __init__(self):
        pass

    def fold(self):
        import alphafold-colabfold
