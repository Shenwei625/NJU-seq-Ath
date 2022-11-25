## 1 About perl 
### 1.1 POD 文档
一种简单而易用的标记型语言
```
POD 文档以=head1开始，=cut结束，=head1前与=cut后添加一空行；perl会忽略POD中的文档
```
+ eg
```perl

=head1 NAME
fetch_fasta.pl -- Get all sequences with the same searched string in a FastA file.
=head1 SYNOPSIS
    perl fetch_fasta.pl -s protein_coding [options]
        Options:
            --help\-h   Brief help message
            --string\-s The sequences we want to fetch
            --file\-f   The FastA file with path
            --stdin     Get FastA from STDIN. It will not been not valid with a provided '--file'
            --rna2dna   Change "U" to "T". Default: False
=cut

```

### 1.2 Perl modular Getopt::Long
用于解析命令行参数的Perl模块

+ GetOptions
```perl
Getopt::Long::GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'string|s=s' => \my $char,
    'file|f=s'   => \my $in_fa,
    'stdin'      => \my $stdin,
    'rna2dna'    => \my $rna2dna,
) or Getopt::Long::HelpMessage(1);

# 'help|h' 接有 | 的选项表示可以简写
# 'string|s=s' 接有 = 的字符串要求接字符串（s）、整数（i）或者浮点（f）等类型的变量
```

### 1.3 Perl modular PerlIO::gzip
读入压缩文件
```perl
open( my $in_fh1, "<:gzip", $ARGV[0] );

#"<:gzip"
```

### 1.4 Shift 函数
这个函数把数组的第一个值移出并且返回它，然后把数组长度减一并且把所有的东西都顺移。如果在数组中不再存在元素，它返回 undef。如果省略了 ARRAY，那么该函数在子过程和格式的词法范围里移动 @_
+ eg
```perl
@arr = ( 1,2,3,4 );
my $str = shift @arr;
# 现在 $str是 1, @arr 是 ( 2,3,4 )
```
+ 变量赋值shitf，接受输入数组的第一个值,然后返回并删除该数组中的第一个项目，一般用于传递给子程序参数

### 1.5 Reverse 函数
改变列表中元素的顺序，并以相反的顺序返回列表
+ eg
```perl
$string = "Hellow World";
print scalar reverse("$string"),"\n";
```
+ 输出
```
dlorW olleH
```

### 1.5 Split 函数
把字符串进行分割并把分割后的结果放入数组中
+ eg
```perl
my @seq1     = split( //, $seq1 );

split /Pattern/, Expression
# Pattern为一个正则表达式，指定了拆分的条件。这里的//与undef用法相同表示将字符串拆分成为每一个字符
```
join 函数与 split 相反，可以把数组组合成字符串
```perl
my $out_seq1 = join( "", @seq1 );
```

### 1.6 cmp 操作符
如果相等,返回0;如果第一个大,返回1;如果第二个大,返回-1

### 1.7 数组前面的$#
他表示数组上最后一个索引的值，如果数组为空，则$#数组为-1

### 1.8 Sort 函数
它对LIST进行排序，并返回排序后的列表，如果指定了SUBNAME，它实际上是个子函数的名字，该子函数对比2个列表元素，并返回一个小于，等于，或大于0的整数，这依赖于元素以何种顺序来sort（升序，恒等，或降序）。也可提供一个BLOCK作为匿名子函数来代替SUBNAME，效果是一样的。
+ eg
```perl
#用法：
#1.sort LIST
#2.sort BLOCK LIST
#3.sort SUBNAME LIST

sort { $map{$a} <=> $map{$b} } keys %map
# 以数字顺序sort

sort { $a cmp $b } @languages
# 以ASCII顺序
```

### 1.9 Perl modular POSIX
统计函数
+ eg
```perl
use POSIX;

POSIX::ceil(3.14)     # 4  向上取整

POSIX::floor(3.14)     # 3   向下取整， 等同于 int(3.14)
```



## 2 About shell
### 2.1 ln 命令
```
ln：link,为某一个文件在另外一个位置建立一个同步的链接。 一种是hard link，又称为硬链接；另一种是symbolic link，又称为符号链接。

-s：对源文件建立符号链接，而非硬链接
-f：强制创建链接，即使目标文件已经存在
```

### 2.2 Time 命令
量测特定指令执行时所需消耗的时间及系统资源等资讯
```bash
time [options] command [arguments...]

# real <=== 实际使用时间（real time）
# user <=== 用户态使用时间（the process spent in user mode）
# sys <=== 内核态使用时间（the process spent in kernel mode）
```

### 2.3 Tee 命令
读取标准输入的数据，并将其内容输出成文件
+ eg
```bash
command [arguments...] | tee test.log
```

### 2.4 awk 命令
与sed一样，均是一行一行的读取，处理；sed作用于一整行的处理，而awk将一行分成数个字段来处理
+ eg
```bash
# 输入：E00516:621:H7KNHCCX2:2:1101:5375:3173   77      *       0       0       *       *       0       0       ATAGTCGATGGACAACAGGTTAGCG       AAFFFJJJJJJJJJJJJJJJJJJJJ       YT:Z:UP
awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'

# 第六列CIGAR有信息（$6!="*"）并且（&&）序列第二次比对时与参开序列一摸一样（$7=="="）
```
+ parallel 命令中的–pipe参数将cat输出分成多个块分派给awk调用,这里每个块处理1千万行数据
```bash
time pigz -dcf output/"${PREFIX}"/rrna.raw.sam.gz |
  parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |'
```


## 3 文件格式
### 3.1 Sam 文件简介
+ eg
```
@HD     VN:1.5  SO:unsorted     GO:query
@SQ     SN:18s  LN:1869
@SQ     SN:28s  LN:5035
@SQ     SN:5-8s LN:157
@PG     ID:bowtie2      PN:bowtie2      VN:2.4.5        CL:"/home/linuxbrew/.linuxbrew/bin/../Cellar/bowtie2/2.4.5/bin/bowtie2-align-s --wrapper basic-0 -p 16 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/hsa_rrna -S output/HeLa_RF_NC/rrna.raw.sam -1 data/HeLa_RF_NC/R1.fq.gz -2 data/HeLa_RF_NC/R2.fq.gz"
E00517:615:HCJYHCCX2:3:1101:6877:1555   99      18s     1123    255     1X20=   =       1123    21      NCAAGGCTGAAACTTAAAGGA   #AAFFJJJJJJJJJJJJJJJJ    AS:i:0  XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:0G20       YS:i:0  YT:Z:CP
E00517:615:HCJYHCCX2:3:1101:6877:1555   147     18s     1123    255     4=1X16= =       1123    -21     GCAANGCTGAAACTTAAAGGA   JJJJ#JJJFJJJJJJJFFFAA    AS:i:0  XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:4G16       YS:i:0  YT:Z:CP
E00517:615:HCJYHCCX2:3:1101:7344:1555   77      *       0       0       *       *       0       0       NGCGGTGAAATGCGTAG       #AAFFJJJJJJJJJJJJ        YT:Z:UP
E00517:615:HCJYHCCX2:3:1101:7344:1555   141     *       0       0       *       *       0       0       CTACGCATTTCACCGCN       AAFFFJJJJJJJJJJJ#        YT:Z:UP

头部区：以'@'开始，体现了比对的一些总体信息。比如比对的SAM格式版本，比对的参考序列，比对使用的软件等。
主体区：比对结果，每一个比对结果是一行，有11个主列和一个可选列。
第一列：QNAME，比对的序列名称，就是fq文件中的read ID，是一条测序read的名称。
第二列：FLAG，比对上的情况
第三列：染色体名称
第四列：POS，比对上的最左边的定位
第五列：MAPQ，比对的质量值。
第六列：CIGAR Extended CIGAR string，M表示匹配、I表示插入、D表示删除、N表示内含子和D类似、S表示替换、H表示剪切、X表示错配。87M表示87个碱基在比对时完全匹配。
第七列：MRNM，这条reads第二次比对的位置，是比对上的参考序列名 。=表示参考序列与reads一模一样，*表示没有完全一模一样的参考序列。
第八列：MPOS，与该reads对应的mate pair reads的比对位置（即mate），若无mate,则为0。
第九列：ISIZE 插入片段长度 例如：200。如果R1端的read和R2端的read能够mapping到同一条Reference序列上（即第三列RNAME相同），则该列的值表示第8列减去第4列加上第6列的值，
第十列：SEQ，和参考序列在同一个链上比对的序列。
第十一列：比对序列的质量（ASCII-33=Phred base quality）reads碱基质量值
```

### 3.2 Fastq 文件格式信息（以二代测序为例）
4行为1个单位
第一行是@开头的解释信息，是这一条read的名字，这个字符串是根据测序时的状态信息转换过来的，中间不会有空格，它是每一条read的唯一标识符，同一份FASTQ文件中不会重复出现，甚至不同的FASTQ文件里也不会有重复；
第二行是序列reads
第三行+开头的信息，在旧版的FASTQ文件中会直接重复第一行的信息，但现在一般什么也不加（节省存储空间）
第四行是对应的测序质量ASCII码
+ eg
```
@E00516:621:H7KNHCCX2:2:1101:1864:3173 2:N:0:ATGTCAAT+AAAGATAA
NCCGCCACCTTCCC
+
#AAFFJJFFJJJJJ
```

### 3.3 Bowtie2 summary文件解读
+ eg
```
50581801 reads; of these:
  50581801 (100.00%) were paired; of these:
    42170200 (83.37%) aligned concordantly 0 times
    7400479 (14.63%) aligned concordantly exactly 1 time
    1011122 (2.00%) aligned concordantly >1 times
    ----
    42170200 pairs aligned concordantly 0 times; of these:
      6058 (0.01%) aligned discordantly 1 time
    ----
    42164142 pairs aligned 0 times concordantly or discordantly; of these:
      84328284 mates make up the pairs; of these:
        84137072 (99.77%) aligned 0 times
        181683 (0.22%) aligned exactly 1 time
        9529 (0.01%) aligned >1 times
16.83% overall alignment rate

不合理比对：
1.比对方向不对，pair-end测序的方向是固定的；
2.read1 和 read2 的插入片段的长度是有限的


结果分为三部分：
第一部分：pair-end模式下，合理比对的结果，read1 和 read2 同时合理的比对到了基因组上
第二部分：pair-end模式下，不合理比对的结果
第三部分：既不能 concordantly,也不能 disconcordantly 1 time 的单端模式比对
```

## 4 一些工具
### 4.1 Bowtie2 简介
序列别对工具，默认使用全局比对，reads 是完整的，看能比对上哪些reference
```
  Read:      GACTGGGCGATCTCGACTTCG
             |||||  |||||||||| |||
  Reference: GACTG--CGATCTCGACATCG
```
+ eg
```bash
bowtie2 -p "${THREAD}" -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/hsa_rrna \
  -1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz \
  -S output/"${PREFIX}"/rrna.raw.sam

# -p:number of alignment threads to launch
# -a:report all alignments; very slow, 输出所有比对结果
# -t:print wall-clock time taken by search phases
# --end-to-end:entire read must align; no clipping (on)
# -D:give up extending after <int> failed extends in a row (15),比对时,将一个种子延长后得到比对结果,如果不产生更好的或次好的比对结果,则该次比对失败.当失败次数连续达到次后,则该条read比对结束. Bowtie2才会继续进行下去
# -R:for reads w/ repetitive seeds, try <int> sets of seeds (2),如果一个read所生成的种子在参考序列上匹配位点过多.当每个种子平均匹配超过300个位置,则通过一个不同的偏移来重新生成种子进行比对. 则是重新生成种子的次数
# -N:max # mismatches in seed alignment; can be 0 or 1 (0)
# -L:length of seed substrings; must be >3, <32 (22),设定种子的长度
# -i:interval between seed substrings w/r/t read len(S,1,1.15),设定两个相邻种子间所间距的碱基数,在--end-to-end模式中默认值为”-iS,1,1.15”.即表示f(x) = 1 + 1.15 *sqrt(x).如果read长度为100,则相邻种子的间距为12.
# --np:penalty for non-A/C/G/Ts in read/ref (1),当匹配位点中read,reference上有不确定碱基(比如N)时所设定的罚分值.Default: 1.
# --xeq:Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record.
# -x:由bowtie2-build所生成的索引文件的前缀
```

### 4.2 cutadapt 简介
去除接头工具
+ eg
```bash
cutadapt -O 6 -m 10 -e 0.1 --discard-untrimmed -j 16 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-o data/${PREFIX}/R1.fq.gz \
-p data/${PREFIX}/R2.fq.gz \
PATH_TO_RAW_DATA/${PREFIX}_R1.fq.gz \
PATH_TO_RAW_DATA/${PREFIX}_R2.fq.gz

# -O:MINLENGTH, --overlap MINLENGTH;如果两端的序列与接头有MINLENGTH个碱基的匹配将会被剔除
# -m:LEN, --minimum-length LEN;如果剔除接头后read长度低于LEN，这条read将会被丢弃
# -e:Maximum allowed error rate
# --discard-untrimmed:丢弃没有接头的read
# -a:3'adapter to be removed from R1
# -A:3'adapter to be removed from R2
# -o:输出 R1 去掉接头的结果
# -p:输出 R2 去掉接头的结果
```

