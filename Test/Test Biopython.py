from Bio import SeqIO
import vcf
import numpy as np

text = np.loadtxt(fname = "binaryTree10.txt")

print(text[2])



# vcf_reader = vcf.Reader(filename="MetaDataFrame.vcf")
# vcf_writer1 = vcf.Writer(open("/Users/leonorschubert/Desktop/modifiedMetaData1.vcf", 'w'), vcf_reader)
#
# record = vcf.model._Record(CHROM=18, POS=5, ID='.', REF=vcf.model._Substitution("G"),
#             ALT=[vcf.model._Substitution("A")], QUAL='.', FILTER='PASS', INFO={},
#             FORMAT=".", sample_indexes=[], samples=None)
# vcf_writer1.write_record(record)
#
#
# vcf_writer1.flush()
#
# vcf_reader2 = vcf.Reader(filename="/Users/leonorschubert/Desktop/modifiedMetaData1.vcf")
# for record2 in vcf_reader2:
#     print (record2.REF)
#
# record = vcf.model._Record(CHROM=18, POS=6, ID='.', REF=vcf.model._Substitution("G"),
#             ALT=[vcf.model._Substitution("T")], QUAL='.', FILTER='PASS', INFO={},
#             FORMAT=".", sample_indexes=[], samples=None)
# vcf_writer1.write_record(record)
# vcf_writer1.flush()
#
# print(":/")
#
# vcf_reader2 = vcf.Reader(filename="/Users/leonorschubert/Desktop/modifiedMetaData1.vcf")
# for record2 in vcf_reader2:
#     print (record2.POS)
#
# print(":/")
#
# vcf_reader2 = vcf.Reader(filename="/Users/leonorschubert/Desktop/modifiedMetaData1.vcf")
# vcf_writer3 = vcf.Writer(open("/Users/leonorschubert/Desktop/modifiedMetaData3.vcf", "w"), vcf_reader)
#
# for record3 in vcf_reader2:
#     print (record3.POS)
#     vcf_writer3.write_record(record3)
#     vcf_writer3.flush()
#
# print("/.")
#
#
# record = vcf.model._Record(CHROM=18, POS=7, ID='.', REF=vcf.model._Substitution("M"),
#             ALT=[vcf.model._Substitution("T")], QUAL='.', FILTER='PASS', INFO={},
#             FORMAT=".", sample_indexes=[], samples=None)
#
# vcf_writer3.write_record(record)
# vcf_writer3.flush()
#
# print("./")
# vcf_reader3 = vcf.Reader(filename="/Users/leonorschubert/Desktop/modifiedMetaData3.vcf")
# for record6 in vcf_reader3:
#     print (record6.POS)
