from Objects import Id,Locus

class Genome():
    dir = '/Users/alexanderfarnum/Documents/Data'
    def __init__(self,gname):
        self.dir += '/'+gname
        fname = self.ext2file(dir=self.dir,fext='gtf')
        from gtfparse import read_gtf
        self.dataframe = read_gtf(fname)
        self.fields = self.dataframe.columns

    def load_datafile(self,fext):
        from Bio import SeqIO
        fname = self.ext2file(dir=self.dir,fext=fext)
        match fext:
            case 'faa' | 'fna':
                return SeqIO.parse(fname,format='fasta')

    def ext2file(self,dir,fext):
        from os import listdir,path
        f_list = listdir(dir)
        f_idx = [x for x in range(len(f_list)) if f_list[x].endswith('.'+fext)]
        return path.join(dir,f_list[f_idx[0]])

    def get_faa(self):
        import os
        from Bio import SeqIO
        f_list = os.listdir(self.dir)
        f_idx = [x for x in range(len(f_list)) if f_list[x].endswith('.faa')]
        f = SeqIO.parse(os.path.join(self.dir,f_list[f_idx[0]]),format='fasta')

        len(f)
        for i,record in f:
            print("ID %s" % record.id)
            print("Sequence length %i" % len(record))


#separate Family into json file
class Family():
    myc = ['Cmyc','Nmyc','Lmyc']
    pur = ['Pura','Purb','Purg1','Purg2']

class Gene():
    """
    gca: RefSeq assembly (NCBI curated version of gcf)
    gcf: GenBank assembly
    """
    def __init__(self,
                 genome):
        self.genome = Genome(gname=genome)

    def get_dna(self,key_type):
        from Objects import Dna
        self.dna = Dna()
        if key_type == 'locus':
            f = self.genome.load_datafile('fna')
            self.locus2entry(file=f,
                             feat=self.dna,
                             locus=self.locus)
        else:           #TODO: add functionality for alt key support
            match key_type:
                case 'gcf':
                    key = self.id.gcf
                case 'hgnc':
                    key = self.id.hgnc
                case _:
                    key_type = 'symbol'
                    key = self.symbol
            self.key2entry(key_type,key)

    def get_rna(self):
        pass

    def get_protein(self,key_type):
        from Objects import Protein
        self.protein = Protein()
        if key_type == 'locus':
            f = super().load_datafile('faa')
            super().locus2entry(file=f,
                                feat=self.protein,
                                locus=self.locus)


    def key2entry(self,key_type,key,feat):
        """
        HGNC test:      dstr = '37102'
        GeneID test:    dstr = '728642'
        """

        delim_map = {
            'ccds':'CCDS:CCDS',
            'geneid':'GeneID:',
            'genbank':'GenBank:',
            'hgnc':'HGNC:HGNC:',
            'mim':'MIM:',
        }

        get_id = lambda key,dstr :[i for i,x in enumerate(self.dataframe[:,self.fields.index('db_xref')]) if split_csv(x)[get_idx(split_csv(x),key)].replace(delim_map[key.lower()],'')==dstr]
        get_idx = lambda lst,key : [i for i,x in enumerate(lst) if x.lower().startswith(key.lower())][0]
        split_csv = lambda s : s.split(',')

        return [self.dataframe[x,:] for x in get_id(key,key_type) if self.dataframe[x,self.fields.index('feature')]==feat][0]
        
    def locus2entry(self,file,feat,locus):
        get_chr = lambda id : True if id == record.id[0:len(id)] else False

        for record in file:
            if get_chr(feat.template+str(locus.chr)):
                print(record.description)
                if record.description.endswith('Primary Assembly'):
                    feat.addentry(seq=record.seq,loc=locus)


class Pura(Gene):
    symbol = 'PURA'
    id = Id(gcf='000001405',hgnc='5813')
    def __init__(self,genome='GRCh37'):
        super().__init__(genome=genome)
        match genome:
            case 'GRCh37':
                self.locus = Locus(chr=5,
                                   start=139493694,
                                   stop=139505204)
            case 'GRCh38':
                self.locus = Locus(chr=5,
                                   start=140114109,
                                   stop=140125619)

class Purb(Gene):
    symbol = 'PURB'
    id = Id(gcf='000001405',hgnc='5813')
    def __init__(self,genome='GRCh37'):
        super().__init__(genome=genome)
        match genome:
            case 'GRCh37':
                self.locus = Locus(chr=5,
                                   start=139493694,
                                   stop=139505204)
            case 'GRCh38':
                self.locus = Locus(chr=5,
                                   start=140114109,
                                   stop=140125619)

class Purg(Gene):
    symbol = 'PURG'
    id = Id(gcf='000001405',hgnc='5813')
    def __init__(self,genome='GRCh37'):
        super().__init__(genome=genome)
        match genome:
            case 'GRCh37':
                self.locus = Locus(chr=5,
                                   start=139493694,
                                   stop=139505204)
            case 'GRCh38':
                self.locus = Locus(chr=5,
                                   start=140114109,
                                   stop=140125619)

class Malinc1(Gene):
    symbol = 'MALINC1'
    id = Id(gcf='000001405',hgnc='5813')
    def __init__(self,genome='GRCh37'):
        super().__init__(genome=genome)
        match genome:
            case 'GRCh37':
                self.locus = Locus(chr=5,
                                   start=139493694,
                                   stop=139505204)
            case 'GRCh38':
                self.locus = Locus(chr=5,
                                   start=140114109,
                                   stop=140125619)