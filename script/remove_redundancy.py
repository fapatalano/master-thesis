from contextlib import ExitStack
import tempfile
from sultan.api import Sultan

def reduce_redundancy(fasta_file):
    with Sultan.load(sudo=False) as s:
        tmp = tempfile.TemporaryDirectory(dir='/tmp')
        clu = tempfile.NamedTemporaryFile(suffix='_cluster')
        DB = tempfile.NamedTemporaryFile(suffix='_DB',dir='/tmp')
        DB_clu_rep = tempfile.NamedTemporaryFile(suffix='_DB_clu_rep',dir='/tmp')
        fasta_DB_clu_rep= tempfile.NamedTemporaryFile(suffix='_fastarepr',dir='/tmp')
        try:
            s.mmseqs(f'createdb {fasta_file} {DB.name}').run()
            s.mmseqs(f'cluster {DB.name} {clu.name} {tmp.name} --min-seq-id 0.3 -c 0.5 --cov-mode 0').run()
            s.mmseqs(f'createsubdb {clu.name} {DB.name} {DB_clu_rep.name}').run()
            s.mmseqs(f'convert2fasta {DB_clu_rep.name} {fasta_DB_clu_rep.name}').run()
        finally:
            clu.close()
            DB.close()
            DB_clu_rep.close()
            fasta_DB_clu_rep.close()


if __name__ == '__main__':
    s = Sultan()
    fasta_file='output.fasta'
    reduce_redundancy(fasta_file)