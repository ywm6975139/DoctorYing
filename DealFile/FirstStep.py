# -*- coding: utf-8 -*-
import os
import sys
#import pandas as pd
import json
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein


reload(sys)
sys.setdefaultencoding('utf8')


#转换文件格式从fastq转换成fasta
def fileConvert(file_source_path,file_target_path):
    records = SeqIO.parse(file_source_path, "fastq")
    count = SeqIO.write(records, file_target_path, "fasta")
    if count>0:
        print "file convert Success! Converted %i records" % count
    else:
        print "file convert Error!"


# 函数功能：输入两段长度为150的基因序列，进行首尾拼接，返回拼接结果。
# 如果没有重复片段，返回 (0,None)
# 如果重复片段小于30，返回 (1,'重复片段')
# 如果重复片段大于等于30，返回 (2,'重复片段')
def seqConcatDeal(seq1, seq2):
    if seq1 == seq2:
        return 2,seq1
    else:
        seqCursor = 1
        while seqCursor < 150:
            if seq1[seqCursor:] == seq2[0:-seqCursor]:
                if seqCursor <= 120:
                    return 2,seq1[:seqCursor] + seq2
                else:
                    return 1,seq1[:seqCursor] + seq2
            else:
                seqCursor = seqCursor + 1
        return 0,None






if __name__ == '__main__':
    file_source_path1="../File/test1.fastq"
    file_target_path1="../File/test1.fasta"
    file_source_path2="../File/test2.fastq"
    file_target_path2="../File/test2.fasta"

    # fileConvert(file_source_path1,file_target_path1)
    # fileConvert(file_source_path2,file_target_path2)


    all_record=[]
    for seq_record in SeqIO.parse(file_target_path1, "fasta"):
        print str(seq_record.id)
        print str(seq_record.seq)
        print str(seq_record.seq.reverse_complement())
        print len(seq_record)


