# -*- coding: utf-8 -*-
import sys
import os
import time
import datetime
import shutil

# import pandas as pd
from Bio import SeqIO

reload(sys)
sys.setdefaultencoding('utf8')


# 第一步：测序数据格式转换，并进行基因片段拼接。
def first_step(file_source_path_1, file_source_path_2, file_target_dir, mk_target_dir_if_not_exist=True):
    print "first_step begin, current time:%s, file_source_path_1:%s, file_source_path_2:%s, file_target_dir:%s" % \
          (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())), file_source_path_1, file_source_path_2,
           file_target_dir)
    start_time = datetime.datetime.now()
    # 校验入参
    if not os.path.exists(file_target_dir):
        if not mk_target_dir_if_not_exist:
            print "file_target_dir is not exist"
            return
        else:
            os.makedirs(file_target_dir)
    if not os.access(file_source_path_1, os.R_OK):
        print "file_source_path_1 is not exist or unable to read"
        return
    if not os.access(file_source_path_2, os.R_OK):
        print "file_source_path_2 is not exist or unable to read"
        return
    # 将fastq格式文件转为fasta格式
    temp_dir = file_target_dir + os.sep + ".temp0"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.mkdir(temp_dir)
    file_target_path_1 = temp_dir + os.sep + os.path.basename(file_source_path_1)[:-1] + "a"
    file_target_path_2 = temp_dir + os.sep + os.path.basename(file_source_path_2)[:-1] + "a"
    if not file_convert(file_source_path_1, file_target_path_1):
        return
    if not file_convert(file_source_path_2, file_target_path_2):
        return
    # 处理序列拼接，生成拼接结果文件，每个文件的序列数不超过50w
    file_path_list = seq_concat_deal(file_target_path_1, file_target_path_2, temp_dir)
    # 合并拼接结果，每个文件的序列数不超过500w
    result_dict = merge_concat_result(file_path_list, file_target_dir)
    print "first_step end,cost time:%ds" % (datetime.datetime.now() - start_time).seconds
    # 删除临时文件
    shutil.rmtree(temp_dir)
    return result_dict


# 转换文件格式从fastq转换成fasta
def file_convert(file_source_path, file_target_path):
    print "file_convert begin, file_source_path:%s, file_target_path:%s" % (file_source_path, file_target_path)
    start_time = datetime.datetime.now()
    records = SeqIO.parse(file_source_path, "fastq")
    count = SeqIO.write(records, file_target_path, "fasta")
    if count > 0:
        print "file convert Success, Converted %i records, cost time:%ds" % \
              (count, (datetime.datetime.now() - start_time).seconds)
        return True
    else:
        print "file convert Error!"
        return False


# 处理基因拼接，生成拼接结果文件，每个文件的序列数不超过50w
def seq_concat_deal(file_target_path_1, file_target_path_2, temp_dir):
    print "seq_concat_deal begin"
    start_time = datetime.datetime.now()
    seq_list = [[], [], []]
    prefix_list = [1, 1, 1]
    suffix_list = ["_none.fasta", "_lt.fasta", "_ge.fasta"]
    file_dict = SeqIO.index(file_target_path_2, "fasta")
    file_path_list = [[], [], []]
    for temp_seq in SeqIO.parse(file_target_path_1, "fasta"):
        seq1 = temp_seq.seq
        seq2 = file_dict[temp_seq.id].seq.reverse_complement()
        repeat_case, concat_result = seq_concat(seq1, seq2)
        temp_seq.seq = concat_result
        seq_list[repeat_case].append(temp_seq)
        if len(seq_list[repeat_case]) >= 50000:
            file_path = temp_dir + os.sep + str(prefix_list[repeat_case]) + suffix_list[repeat_case]
            print "create split file, path:%s" % file_path
            SeqIO.write(seq_list[repeat_case], file_path, "fasta")
            file_path_list[repeat_case].append(file_path)
            seq_list[repeat_case] = []
            prefix_list[repeat_case] += 1
    for index in range(2):
        if len(seq_list[index]) > 0:
            file_path = temp_dir + os.sep + str(prefix_list[index]) + suffix_list[index]
            print "create split file, path:%s" % file_path
            SeqIO.write(seq_list[index], file_path, "fasta")
            file_path_list[index].append(file_path)
    print "seq_concat_deal end, cost time:%ds" % (datetime.datetime.now() - start_time).seconds
    return file_path_list


# 输入两段基因序列，进行首尾拼接，返回拼接结果。
def seq_concat(seq1, seq2):
    bytes1 = bytes(seq1)
    bytes2 = bytes(seq2)
    len1 = len(bytes1)
    len2 = len(bytes2)
    index = len1 - len2 if len1 > len2 else 0
    while index < len1:
        if bytes1[index] == bytes2[0]:
            index2 = index + 1
            match = True
            while index2 < len1:
                if bytes1[index2] != bytes2[index2 - index]:
                    match = False
                    break
                index2 += 1
            if match is True:
                match_count = len1 - index
                if match_count >= 30:
                    return 2, seq1 + seq2[match_count:]
                return 1, seq1 + seq2[match_count:]
        index += 1
    return 0, seq1 + seq2


# 合并拼接结果，每个文件的序列数不超过500w
def merge_concat_result(file_path_list, file_target_dir):
    print "merge_concat_result begin"
    start_time = datetime.datetime.now()
    suffix_list = ["_none.fasta", "_lt.fasta", "_ge.fasta"]
    return_list = [[], [], []]
    return_dict = {'none': return_list[0], 'lt': return_list[1], 'ge': return_list[2]}
    for index in range(3):
        path_list = file_path_list[index]
        file_count = len(path_list)
        if file_count == 0:
            continue
        prefix = 1
        count = 0
        file_name = str(prefix) + suffix_list[index]
        file_path = file_target_dir + os.sep + file_name
        print "create result file, path:%s" % file_path
        result_file = open(file_path, 'w')
        for split_file_path in path_list:
            for line in open(split_file_path):
                result_file.writelines(line)
            count += 1
            if count % 10 == 0 and count < file_count:
                result_file.close()
                return_list[index].append(file_name)
                prefix += 1
                file_name = str(prefix) + suffix_list[index]
                file_path = file_target_dir + os.sep + file_name
                print "create result file, path:%s" % file_path
                result_file = open(file_path, 'w')
        result_file.close()
        return_list[index].append(file_name)
    print "merge_concat_result end, cost time:%ds, return_dict:%s" % \
          ((datetime.datetime.now() - start_time).seconds, return_dict)
    return return_dict


# 构造测试数据
def make_test_file(path1, path2, count=10000000):
    start_time = datetime.datetime.now()
    line3 = "+\n"
    line4 = "#AAFFJJJJJJJJJJJ7JJJJJJJJ<JJ-AJJJJJJJJFF7FJFAJJJJJJJJJJFJJ7AF#<7AAJ<FFJJJJAJAJAAFAAF#<FJJJFF#7<F)-77JJ<7J7)<<AAJFJAF<<A<))7<F<77)A)7AFJJJ-7-----7))7)<\n"
    seq11 = "NTGTCCCTCACCTGCGCTGTCTATGGTGGGTCCTTAAGTGGCTACTACTGGAGCTGGGTCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGGCGAGTC\n"
    seq12 = "NGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCGTAGTAGTGGTGGGACCAGCTGCTGNTATACCCTTTTCCCGCACAGTANTACACGGNCGTGTCCGCAGCGGTCACAGAGCTCAGCGTCAGGGAGAACTGGAACTTGGACGTGGC\n"
    seq21 = "NGTGACCATTGTCCCTTGGCCCCACATATCAAAAGCATCATCGCCGATAACTCCCCCATAGGTAATCGGGCAGTGTGCACAATAATATGTGCCTGTGTCCTCAGGGTCCATGTTGGTCATTGTTATGACCACCTGGTTTTTGGAGGTGTC\n"
    seq22 = "NTCACGCTGACCTGCACCTTCTCTGGCTTCTCACTNACCACTCATGGAGTGGGTGTGGGCTNGATCCGNCAGCCCCCAGGAAAGNCCCTGGANTGGCTTGCACNCATTTATTGGGNTAGTGATAAGCGCTACAGCCCATCTCTGAAGAAC\n"
    seq31 = "NTGTCCCTCACCTGCACTGTCTCTGGTGACTCCATCAGTAGTTATTACTGGAGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCGCCCTCAAGAGTCGAGGT\n"
    seq32 = "NGTGACCATTGTCCCTTGGCCCCAGATAAAAGGGGNGGTCAGCGTAGTTCGGCTCGCACAGNAATACANGGCCGTGTCCGCAGCNGTCACAGNGCTCAGCTTCNGGGAGATCTGGNTCTTGGACGTGTCTACGGATATCGTAACTCGACT\n"
    seq41 = "NTCACGCTGACCTGCACCTTCTCTGGCTTCTCACTCACCACTCGTGGAGTGGGTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGGCCCTGGAGTGGCTTGCACTCATTTATTGGAATGATGATAAGCGCTACAGCCCATCTCTGAAGAGC\n"
    seq42 = "NGTGACCATTGTCCCTTGGCCCCAGACAAATTTCCNCCCACGTCTGTGTGCACAGTAATATNTGGCTGNGTCCACAGGGTCCATNTTGGTCANTGTAANGACCNCCTGGTTTTTGNAGGTGTCCTTGGTGATGGTGAGCCTGCTCTTCAG\n"
    seq51 = "NTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATTAATAGTTACTATTGGATCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCGTATATCTATTACAGTGGGATCACCAACTACAACCCCGCCCTCAAGAGCCTAGTC\n"
    seq52 = "NGTGACCGTGGTCCCTTGGCCCCAGATATCAAAAGCATCCCTTACGTAGTTCCAGTCCTTCCTCGCACAGTAATACACGGCCGTNTCCGCGGNGGTCACAGAGCTCAGCTTCAGGGAGAACTGGTTCTTGGACGTGGCGACTGAAATGAT\n"
    list1 = [seq11, seq21, seq31, seq41, seq51]
    list2 = [seq12, seq22, seq32, seq42, seq52]
    file1 = open(path1, 'w')
    file2 = open(path2, 'w')
    range_count = count / 5
    for index in range(range_count):
        fill = str(index).zfill(7)
        prefix = "@J00103:57:H53CJBBXX:5:11" + fill[0:2] + ":4" + fill[2:5] + ":1" + fill[5:]
        id1 = prefix + "1 1:N:0:NGATGT\n"
        id2 = prefix + "2 1:N:0:NGATGT\n"
        id3 = prefix + "3 1:N:0:NGATGT\n"
        id4 = prefix + "4 1:N:0:NGATGT\n"
        id5 = prefix + "5 1:N:0:NGATGT\n"
        ids = [id1, id2, id3, id4, id5]
        for x in range(5):
            file1.writelines(ids[x])
            file1.writelines(list1[x])
            file1.writelines(line3)
            file1.writelines(line4)
            file2.writelines(ids[x])
            file2.writelines(list2[x])
            file2.writelines(line3)
            file2.writelines(line4)
    file1.close()
    file2.close()
    print "make_test_file end,cost time:%ds" % (datetime.datetime.now() - start_time).seconds


if __name__ == '__main__':
    # make_test_file("D:\\work\\DoctorYingTest\\1000w1.fastq", "D:\\work\\DoctorYingTest\\1000w2.fastq", 10000000)
    first_step("D:\\work\\DoctorYingTest\\1w1.fastq", "D:\\work\\DoctorYingTest\\1w2.fastq",
               "D:\\work\\DoctorYingTest\\result")
