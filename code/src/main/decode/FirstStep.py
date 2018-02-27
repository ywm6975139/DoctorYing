# -*- coding: utf-8 -*-
import sys
import os
import datetime
import shutil

# import pandas as pd
from Bio import SeqIO

reload(sys)
sys.setdefaultencoding('utf8')


# 第一步：测序数据格式转换，并进行基因片段拼接。
# 入参：两个fastq格式文件路径、输出文件夹路径、输出路径不存在时退出还是创建文件夹
def first_step(file_source_path_1, file_source_path_2, file_target_dir, mk_target_dir_if_not_exist=True):
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
    temp_dir = file_target_dir + os.sep + "temp"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.mkdir(temp_dir)
    file_target_path_1 = temp_dir + os.sep + os.path.basename(file_source_path_1)[:-1] + "a"
    file_target_path_2 = temp_dir + os.sep + os.path.basename(file_source_path_2)[:-1] + "a"
    print "file_source_path_1:%s" % file_source_path_1
    print "file_source_path_2:%s" % file_source_path_2
    print "file_target_path_1:%s" % file_target_path_1
    print "file_target_path_2:%s" % file_target_path_2
    print "file_target_dir:%s" % file_target_dir
    if not file_convert(file_source_path_1, file_target_path_1):
        return
    if not file_convert(file_source_path_2, file_target_path_2):
        return

    # 将第一个fasta文件进行分割，每个文件不超过100,000条序列
    file_count = split_target_file(file_target_path_1, temp_dir)
    print "file_count:%s" % file_count

    # 处理序列拼接
    for index in range(1, file_count + 1):
        seq_concat_deal(index, temp_dir, file_target_path_2)

    # 合并拼接结果文件，再分割用于上传到网站
    merge_split_result_file(temp_dir, file_count)
    none_file_path_list = split_upload_file(temp_dir + os.sep + "none.fasta", file_target_dir, "_none.fasta")
    lt_file_path_list = split_upload_file(temp_dir + os.sep + "lt.fasta", file_target_dir, "_lt.fasta")
    ge_file_path_list = split_upload_file(temp_dir + os.sep + "ge.fasta", file_target_dir, "_ge.fasta")
    return_dict = {'none': none_file_path_list, 'lt': lt_file_path_list, 'ge': ge_file_path_list}
    shutil.rmtree(temp_dir)
    print "first_step end"
    return return_dict


# 转换文件格式从fastq转换成fasta
def file_convert(file_source_path, file_target_path):
    print "file_convert begin, file_source_path:%s, file_target_path:%s" % (file_source_path, file_target_path)
    records = SeqIO.parse(file_source_path, "fastq")
    count = SeqIO.write(records, file_target_path, "fasta")
    if count > 0:
        print "file convert Success! Converted %i records" % count
        return True
    else:
        print "file convert Error!"
        return False


# 将第一个fasta文件进行分割，每个文件不超过100,000条序列，返回分割所得文件数
def split_target_file(file_target_path, temp_dir):
    print "split_target_file begin"
    count = 1
    current_size = 0
    split_file_path = temp_dir + os.sep + str(count) + "_src.fasta"
    file_list = []
    for temp_seq in SeqIO.parse(file_target_path, "fasta"):
        file_list.append(temp_seq)
        current_size = current_size + 1
        if current_size >= 100000:
            print "create src split file,path:%s" % split_file_path
            SeqIO.write(file_list, split_file_path, "fasta")
            file_list = []
            current_size = 0
            count = count + 1
            split_file_path = temp_dir + os.sep + str(count) + "_src.fasta"
    if len(file_list) > 0:
        print "create src split file,path:%s" % split_file_path
        SeqIO.write(file_list, split_file_path, "fasta")
    else:
        count = count - 1
    print "split_target_file end, result file count is:", count
    return count


# 处理基因拼接
def seq_concat_deal(index, temp_dir, file_target_path):
    print "seq_concat_thread begin, index:", index
    src_file_path = temp_dir + os.sep + str(index) + "_src.fasta"
    none_file_path = temp_dir + os.sep + str(index) + "_none.fasta"
    lt_file_path = temp_dir + os.sep + str(index) + "_lt.fasta"
    ge_file_path = temp_dir + os.sep + str(index) + "_ge.fasta"
    none_seq_list = []
    lt_seq_list = []
    ge_seq_list = []
    file_dict = SeqIO.index(file_target_path, "fasta")
    for temp_seq in SeqIO.parse(src_file_path, "fasta"):
        seq1 = temp_seq.seq
        seq2 = file_dict[temp_seq.id].seq.reverse_complement()
        repeat_case, concat_result = seq_concat(seq1, seq2)
        temp_seq.seq = concat_result
        if 0 == repeat_case:
            none_seq_list.append(temp_seq)
        elif 1 == repeat_case:
            lt_seq_list.append(temp_seq)
        else:
            ge_seq_list.append(temp_seq)
    file_dict.close()
    if len(none_seq_list) > 0:
        print "create none_file, path:%s" % none_file_path
        SeqIO.write(none_seq_list, none_file_path, "fasta")
    if len(lt_seq_list) > 0:
        print "create lt_file, path:%s" % lt_file_path
        SeqIO.write(lt_seq_list, lt_file_path, "fasta")
    if len(ge_seq_list) > 0:
        print "create ge_file, path:%s" % ge_file_path
        SeqIO.write(ge_seq_list, ge_file_path, "fasta")
    print "seq_concat_thread end, index:", index


# 输入两段基因序列，进行首尾拼接，返回拼接结果。
# 如果没有重复片段，返回 (0,'拼接结果')
# 如果重复片段小于30，返回 (1,'拼接结果')
# 如果重复片段大于等于30，返回 (2,'拼接结果')
def seq_concat(seq1, seq2):
    if seq1 == seq2:
        return 2, seq1
    else:
        len1 = len(seq1)
        len2 = len(seq2)
        seq_cursor = 1
        if len1 > len2:
            if seq1[len1 - len2:] == seq2:
                return 2, seq1
            seq_cursor = len1 - len2 + 1
        elif seq1 == seq2[:len1 - len2]:
            return 2, seq2
        while seq_cursor < len1:
            if seq1[seq_cursor:] == seq2[0:-seq_cursor]:
                if len1 - seq_cursor >= 30:
                    return 2, seq1[:seq_cursor] + seq2
                else:
                    return 1, seq1[:seq_cursor] + seq2
            else:
                seq_cursor = seq_cursor + 1
        return 0, seq1 + seq2


# 合并拼接结果文件
def merge_split_result_file(temp_dir, file_count):
    print "merge_split_result_file begin"
    none_file = open(temp_dir + os.sep + "none.fasta", 'w')
    lt_file = open(temp_dir + os.sep + "lt.fasta", 'w')
    ge_file = open(temp_dir + os.sep + "ge.fasta", 'w')
    for index in range(1, file_count + 1):
        none_split_file_path = temp_dir + os.sep + str(index) + "_none.fasta"
        lt_split_file_path = temp_dir + os.sep + str(index) + "_lt.fasta"
        ge_split_file_path = temp_dir + os.sep + str(index) + "_ge.fasta"
        if os.path.exists(none_split_file_path):
            for line in open(none_split_file_path):
                none_file.writelines(line)
        if os.path.exists(lt_split_file_path):
            for line in open(lt_split_file_path):
                lt_file.writelines(line)
        if os.path.exists(ge_split_file_path):
            for line in open(ge_split_file_path):
                ge_file.writelines(line)
    none_file.close()
    lt_file.close()
    ge_file.close()
    print "merge_split_result_file end"


# 分割待上传到网站文件
def split_upload_file(file_path, file_target_dir, file_name_suffix):
    print "split_upload_file begin, file_path:%s" % file_path
    file_len = len(list(SeqIO.parse(file_path, "fasta")))
    print "file_len:", file_len
    if 0 == file_len:
        return
    count = 1
    current_size = 0
    result_file_path = file_target_dir + os.sep + str(count) + file_name_suffix
    result_file_path_list = []
    seq_list = []
    for temp_seq in SeqIO.parse(file_path, "fasta"):
        seq_list.append(temp_seq)
        current_size = current_size + 1
        if current_size >= 500000:
            print "create split upload file, path:%s" % result_file_path
            SeqIO.write(seq_list, result_file_path, "fasta")
            result_file_path_list.append(str(count) + file_name_suffix)
            seq_list = []
            current_size = 0
            count = count + 1
            result_file_path = file_target_dir + os.sep + str(count) + file_name_suffix
    if len(seq_list) > 0:
        print "create split upload file, path:%s" % result_file_path
        SeqIO.write(seq_list, result_file_path, "fasta")
        result_file_path_list.append(str(count) + file_name_suffix)
    else:
        count = count - 1
    print "split_upload_file end, result file count:", count
    return result_file_path_list


if __name__ == '__main__':
    start_time = datetime.datetime.now()
    first_step("D:\\work\\DoctorYingTest\\test1.fastq", "D:\\work\\DoctorYingTest\\test2.fastq", "D:\\work\\DoctorYingTest\\result")
    end_time = datetime.datetime.now()
    print "cost time:%ds" % (end_time - start_time).seconds
    # file_source_path1 = "../File/test1.fastq"
    # file_target_path1 = "../File/test1.fasta"
    # file_source_path2 = "../File/test2.fastq"
    # file_target_path2 = "../File/test2.fasta"
    # # file_convert(file_source_path1,file_target_path1)
    # # file_convert(file_source_path2,file_target_path2)
    # for seq_record in SeqIO.parse(file_target_path1, "fasta"):
    #     print str(seq_record.id)
    #     print str(seq_record.seq)
    #     print str(seq_record.seq.reverse_complement())
    #     print len(seq_record)
