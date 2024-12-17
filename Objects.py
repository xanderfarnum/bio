

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
    template='NC_00000'
    def __init__(self):
        pass

    def addentry(self,seq,loc):
        self.start = loc.start
        self.stop = loc.stop
        self.seq = seq[self.start:self.stop]


class Rna(Object):
    def __init__(self):
        pass

class Protein(Object):
    template='lcl|NC_00000'
    def __init__(self):
        pass

    def fold(self):
        import alphafold

class Id():
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            setattr(self,key,value)

class Locus():
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            setattr(self,key,value)




if __name__ == "__main__":
    from Genes import *

    gene = Pura()
    gene.get_dna(key_type='locus')

