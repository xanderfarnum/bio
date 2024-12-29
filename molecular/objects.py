class Object():
    def __init__(self):
        pass

    def set_seq_idx(self,start,stop):
        self.start = start
        self.stop = stop

    def get_start(self):
        """
        Finds first start codon within sequence
        """
        return self.seq.upper().find(sub=self.cstart)

    def get_stop(self,start):
        """
        Find first stop codon within sequence. Must move in accordance with triplets
        TODO: Iff premature stop codon, return variant and WT
        """
        import numpy as np
        trip_seq = np.arange(start=start+3,
                             stop=len(self.seq)-start,
                             step=3)
        for i in trip_seq: #range(len(trip_seq)-1): #should encounter stop codon prior to ultimate element of trip_seq
            if self.seq[i:i+3].upper() in self.cstop:
                return i

class Dna(Object):
    template='NC_00000'
    cstart='ATG'
    cstop=['TAA','TAG','TGA']
    def __init__(self,rna=None):
        if rna:
            from Bio.Seq import back_transcribe
            self.seq = back_transcribe(rna.seq)

    def set_seq(self,seq,start,stop):
        super().set_seq_idx(start,stop)
        self.seq = seq[self.start:self.stop].upper()

        from Bio.SeqUtils import gc_fraction
        self.gc = gc_fraction(self.seq)

    def get_seq(loc):
        self = Dna()


class Rna(Object):
    cstart='AUG'
    cstop=['UAA','UAG','UGA']
    def __init__(self,dna=None):
        if dna:
            from Bio.Seq import transcribe
            self.seq = transcribe(dna.seq)

    def set_seq(self,seq,start,stop):
        super().set_seq_idx(start,stop)
        self.seq = seq[self.start:self.stop]

class Protein(Object):
    template='lcl|NC_00000'
    def __init__(self,
                 seq=None):
        if seq:     #Instantiate protein object from existing Dna or Rna sequence
            from Bio.Seq import translate
            self.start = seq.get_start()
            self.stop = seq.get_stop(start=self.start)
            self.seq = translate(seq.seq[self.start:self.stop])

    def set_seq(self,seq,start,stop):
        super().set_seq_idx(start,stop)
        self.seq = seq[self.start:self.stop]        

    def fold(self):
        import alphafold



class Id():
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            setattr(self,key,value)

class Location():
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            setattr(self,key,value)
