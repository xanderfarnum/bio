import os
# from gtfparse.read_gtf import parse_gtf,parse_gtf_pandas

from Bio import SeqIO
from abc import ABC,abstractmethod

class AbsFile(ABC):
    dir = {
    'base': '/Users/alexanderfarnum/Documents',
    'data': 'Data',
    'code': 'Code',
    }

    def __init__(self,file_name):
        self.file_name = file_name

class Fasta(AbsFile):
    def __init__(self,file_name,species='Human'):
        """
        filename.fa: nucleic acid or amino acid sequence
        filename.fna: nucleic acid sequence
        filename.faa: amino acid sequence
        """
        super().__init__(file_name)
        file_name,ext = file_name.split('.')

        match ext:
            case 'fna':
                self.nucletoide = True
            case 'faa':
                self.protein = True

        self.fpath = os.path.join(self.dir['base'],self.dir['data'],'RefSeq',species,self.file_name)

        match ext:
            case 'fa' | 'fna' | 'faa':
                ftype = 'fasta'

        self.fiter = SeqIO.parse(open(self.fpath),ftype)

    def iterate(self,idx:list):
        for record in self.fiter:
            print(record.id[idx])

class Annotations(AbsFile):
    ext = 'gtf'
    def __init__(self,file_name,species='Human',gene=None):
        super().__init__(file_name)
        self.fpath = os.path.join(self.dir['base'],self.dir['data'],'RefSeq',species,self.file_name+'.'+self.ext)

        # self.data = parse_gtf(self.fpath)
        # idx = self._match_gene_name(data,gene)
        # pos = self._get_gene_pos(data,idx)

        # attr_idx = self.data.columns.index("attribute")
        # for i in range(self.data.height):
        #     pos = self._match_gene_name(self.data[1,attr_idx],gene=gene)
        #     if pos is None:
        #         continue
        #     else:
        #         break

    @staticmethod
    def _match_gene_name(data,gene):
        col_idx = data.columns.index("attribute")
        gene_name = 'gene_name'

        for i in range(data.height):
            start_idx = data[i,col_idx].find(gene_name)
            substr = data[i,col_idx][start_idx+len(gene_name):]
            end_idx = substr.find(';')
            substr = substr[:end_idx]
            substr = substr.replace('\"','').replace('\'','').strip()
            if substr == gene.upper():
                t = []
                # return i
        print('gene name not found in annotations file')
        input("Press Enter to continue...")
            
    @staticmethod
    def _get_gene_pos(data,row_idx):
        col_idx = [data.columns.index("start"),data.columns.index("end")]

        

class _Variant(AbsFile):
    def __init__(self,file_dir,file_name):
        super().__init__(file_name)
        self.fpath = os.path.join(self.dir['base'],self.dir['data'],'Variants',file_dir,file_name)

    @abstractmethod
    def read(self):
        pass

    def blast(self,seq,ref):
        pass

class Vcf(_Variant):
    def __init__(self,
                 file_name,
                 gene,
                 species='Human',
                 chrom=None,
                 variant='all',
                 ref_genome=None,
                 annotations=None):
        super().__init__(file_name)
        self.data, self.headers = self.read()

        if species is not None and gene is not None:
            self.gene = {
                'name':gene,
                'sequence':Gene(species,gene)
            }

            if ref_genome is not None:
                if annotations is not None:
                    self.genome = {
                        'name':ref_genome,
                        'sequence':Fasta(ref_genome),
                        'annotations':Annotations(annotations,gene=gene)
                    }
                else:
                    self.genome = {
                        'name':ref_genome,
                        'sequence':Fasta(ref_genome)
                    }


    def read(self):
        num_header = 0
        with open(self.fpath) as f:
            for line in f.readlines():
                if line.startswith("##"):
                    num_header += 1
                else:
                    break
        from pandas import read_csv
        data = read_csv(self.fpath, sep="\t", skiprows=num_header)
        data = data.rename({"#CHROM": "CHROM"}, axis=1)
        headers = list(data.columns.values)
        return data, headers
    
    def variant(self,species='Human',gene=None,chrom=None,variant='all',ref_genome=None,annotations=None):
        if species is not None and gene is not None:
            self.gene = {
                'name':gene,
                'sequence':Gene(species,gene)
            }

            if ref_genome is not None:
                if annotations is not None:
                    self.genome = {
                        'name':ref_genome,
                        'sequence':Fasta(ref_genome),
                        'annotations':Annotations(annotations,gene=gene)
                    }
                else:
                    self.genome = {
                        'name':ref_genome,
                        'sequence':Fasta(ref_genome)
                    }
        t = []



        # if chrom is None and variant is None:
        #     return self.data
        # elif chrom is not None:
        #     idx = []
        #     chrom_set = set(self.data['CHROM'])
        #     for i in chrom_set:
        #         idx.append(list(self.data['CHROM']).index(i))
        #     order = [i for i, x in sorted(enumerate(idx), key=lambda x: x[1])]
        #     self.chrom_idx = [[list(chrom_set)[i] for i in order],[idx[i] for i in order]]
        #     start_idx = [x == 'chr'+str(chrom) for x in self.chrom_idx[0]].index(True)
        #     data_subset = self.data[self.chrom_idx[1][start_idx]:self.chrom_idx[1][start_idx+1]]
        #     if variant is not None:
        #         buffer = 200    #number of base pairs to include on either side of variant
        #         if genome is None:
        #             pass
        #         else:
        #             gfile = Genome(genome)
        #             #get appropriate idx, read it with variant substitution and as original genome
        #     else:   #intra-chromosomal bp position specified
        #         print(data_subset)

        #     return data_subset
        # elif variant is not None:    #absolute genome bp position specified
        #     pass

class Cram(_Variant):
    """
    Compressed Reference-oriented Alignment Map
    """
    def __init__(self,file_dir='C',file_name='3004831',ref='GRCh37'):
        super().__init__(file_dir,file_name)

        match ref:
            case 'GRCh37':
                ref = 'GRCh37/GRCh37.p13'
        import pysam
        self.file = pysam.AlignmentFile(filename=self.fpath+'.cram',
                                        mode="rc",
                                        reference_filename=os.path.join(self.dir['base'],self.dir['data'],ref+'.faa'))
        self.ref = Fasta(ext='fna')
        data = self.read(1,1000,2000)
        self.file.close()


    def read(self,chrom,start,stop):
        for read in self.file.fetch("chr"+str(chrom), start, stop):
            print(read.query_name)


class Bam(_Variant):
    def __init__(self,file_name='test'):
        super().__init__(file_name)
    
    def read(self):
        import pysam
        self.file = pysam.AlignmentFile(filename=os.path.join(self.dir['base'],self.dir['data'],self.file_name+'.bam'),
                                        mode="rb",)
                                        # reference_filename=os.path.join(self.dir['base'],self.dir['data'],ref+'.faa'))


    def write_sample(self):
        import pysam

        header = { 'HD': {'VN': '1.0'},
                    'SQ': [{'LN': 1575, 'SN': 'chr1'},
                        {'LN': 1584, 'SN': 'chr2'}] }

        with pysam.AlignmentFile('tmpfile.bam', "wb", header=header) as outf:
            a = pysam.AlignedSegment()
            a.query_name = "read_28833_29006_6945"
            a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            a.flag = 99
            a.reference_id = 0
            a.reference_start = 32
            a.mapping_quality = 20
            a.cigartuples = ((0,10), (2,1), (0,25))
            a.next_reference_id = 0
            a.next_reference_start=199
            a.template_length=167
            a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
            a.tags = (("NM", 1),
                    ("RG", "L1"))
            outf.write(a)


if __name__ == "__main__":
    # vcf_fname = 'clinvar_20240917.vcf'
    # gene = 'pura'
    # genome = 'hg19'
    # genome_annotations = 'gencode_basic'

    # vcf1 = Vcf(file_name=vcf_fname,
    #            gene=gene,
    #            ref_genome=genome,
    #            annotations=genome_annotations)
    
    t = Cram()
