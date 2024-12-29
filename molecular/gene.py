class Gene():
    """
    gca: RefSeq assembly (NCBI curated version of gcf)
    gcf: GenBank assembly
    """
    def __init__(self,obj,source):
        match source:
            case 'remote':      # load via Entrez
                from Bio import Entrez
                from molecular.utils import Query
                if obj.parameters.species: 
                    q = Query(str=obj.parameters.species+"[Orgn]",delim='AND')
                if obj.parameters.gene: 
                    q.append(str=obj.parameters.gene+"[Gene]")

                stream = Entrez.esearch(db="nucleotide",
                                        term=q.construct(),
                                        idtype="acc")
                self.data = Entrez.read(stream)
            case 'local':
                from molecular.genome import Genome
                from molecular.utils import get_subdata, list2attrs
                from molecular.objects import Dna

                self.genome = Genome(gname=obj.parameters.genome)
                self.name = obj.parameters.gene
                list2attrs(obj=self,_dict=get_subdata(obj))
                self.dna = Dna.get_seq()


    def entry_from_file(self,
                        key_type,
                        feat,
                        fext):
        if key_type == 'position':
            f = self.genome.load_datafile(fext)
            self.genome.position2entry(file=f,
                                       feat=feat,
                                       position=self.position['gene']) #TODO: make dynamic reference if protein calls sme super function
        else:           #TODO: add functionality for alt key support
            match key_type:
                case 'gcf':
                    key = self.id.gcf
                case 'hgnc':
                    key = self.id.hgnc
                case _:
                    key_type = 'symbol'
                    key = self.symbol
            self.genome.key2entry(key_type,key)

    def get_dna(self,key_type=None):
        from molecular.objects import Dna
        if hasattr(self,'rna'):
            self.dna = Dna(rna=self.rna)    #back transcribe from rna
        else:
            self.dna = Dna()
            self.entry_from_file(key_type=key_type,
                                 feat=self.dna,
                                 fext='fna')

    def get_rna(self,key_type=None):
        from molecular.objects import Rna
        if hasattr(self,'dna'):
            self.rna = Rna(dna=self.dna)    #back transcribe from rna
        else:
            self.rna = Rna()
            self.entry_from_file(key_type=key_type,
                                 feat=self.rna,
                                 fext='fna')

    def get_protein(self,key_type=None):
        from molecular.objects import Protein
        if hasattr(self, 'dna'):
            self.protein = Protein(seq=self.dna)    #translate from dna
        elif hasattr(self, 'rna'):
            self.protein = Protein(seq=self.rna)    #translate from rna
        else:   # get from file
            self.protein = Protein()
            if key_type == 'position':
                f = super().load_datafile('faa')
                self.genome.position2entry(file=f,
                                    feat=self.protein,
                                    position=self.position)