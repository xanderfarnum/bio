from Gene import *
from Objects import print_seq

object = Pura()
object.get_dna(key_type='position')
object.get_rna()
object.get_protein()

print_seq(object.protein.seq)
