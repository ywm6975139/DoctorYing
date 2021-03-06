# _*_ coding:utf-8 _*_
# 输入两段基因序列，进行首尾拼接，返回拼接结果。
# 如果没有重复片段，返回 (0,'拼接结果')
# 如果重复片段小于30，返回 (1,'拼接结果')
# 如果重复片段大于等于30，返回 (2,'拼接结果')
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


# 测试两个完全相同的序列
seq1 = 'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
seq2 = 'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
if seq_concat(seq1, seq2) == (2, seq1):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]
# 测试重复片段大于30的两个序列
seq1 = 'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
seq2 = 'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATAMKTESTDATA'
if seq_concat(seq1, seq2) == (2,
                             'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATAMKTESTDATA'):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]
# 测试重复片段小于30的两个序列
seq1 = 'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
seq2 = 'CKTESTDATAMKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
if seq_concat(seq1, seq2) == (1,
                             'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATAMKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]
# 测试重复片段等于30的两个序列
seq1 = 'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
seq2 = 'CKTESTDATACKTESTDATACKTESTDATAMKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
if seq_concat(seq1, seq2) == (2,
                             'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATAMKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]
# 测试没有重复片段的两个序列
seq1 = 'CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
seq2 = 'MKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA'
if seq_concat(seq1, seq2) == (0, "CKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATAMKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATACKTESTDATA"):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]

seq1 = 'ABCDEFGH'
seq2 = 'ABCDEFGH'
if seq_concat(seq1, seq2) == (1, "ABCDEFGH"):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]

seq1 = 'ABCDEF'
seq2 = 'ABCDEFGH'
if seq_concat(seq1, seq2) == (1, "ABCDEFGH"):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]


seq1 = 'ABCDEFGH'
seq2 = 'DEFGH'
if seq_concat(seq1, seq2) == (1, "ABCDEFGH"):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]

seq1 = 'ABCDEF'
seq2 = 'BCDEFGH'
if seq_concat(seq1, seq2) == (1, "ABCDEFGH"):
    print 'pass'
else:
    print 'error'
    print seq_concat(seq1, seq2)[0]
    print seq_concat(seq1, seq2)[1]