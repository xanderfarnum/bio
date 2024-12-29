
def get_from_entrez(db=None):
    from Bio import Entrez
    Entrez.email = "ajfarnum@yahoo.com"
    if not db:
        stream = Entrez.einfo()
        record = Entrez.read(stream)
        result = record["DbList"]
    else:
        stream = Entrez.einfo(db=db)
        record = Entrez.read(stream)
        result = record["DbInfo"]["Description"]+'\n'
        for field in record["DbInfo"]["FieldList"]:
            result+='\t'+("%(Name)s -- %(FullName)s -- %(Description)s"%field)+'\n'
    return result
