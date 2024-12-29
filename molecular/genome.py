class Genome():
    dir = '/Users/alexanderfarnum/Documents/Data'
    def __init__(self,gname):
        self.dir += '/'+gname
        fname = self._ext2file(dir=self.dir,fext='gtf')
        from gtfparse import read_gtf
        self.dataframe = read_gtf(fname)
        self.fields = self.dataframe.columns

    def _ext2file(self,dir,fext):
        from os import listdir,path
        f_list = listdir(dir)
        f_idx = [x for x in range(len(f_list)) if f_list[x].endswith('.'+fext)]
        return path.join(dir,f_list[f_idx[0]])


    def load_datafile(self,fext):
        from Bio import SeqIO
        fname = self.ext2file(dir=self.dir,fext=fext)
        match fext:
            case 'faa' | 'fna':
                return SeqIO.parse(fname,format='fasta')

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
        
    def loc2entry(self,file,feat,loc):
        get_chr = lambda id : True if id == record.id[0:len(id)] else False

        for record in file:
            if get_chr(feat.template+str(loc.chr)):
                print(record.description)
                if record.description.endswith('Primary Assembly'):
                    feat.set_seq(seq=record.seq,
                                 start=loc.start,
                                 stop=loc.stop)