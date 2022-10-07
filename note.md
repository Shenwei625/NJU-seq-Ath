# 2022.9.21
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



# 2022.9.28
### 1 ln命令
```
ln：link,为某一个文件在另外一个位置建立一个同步的链接。 一种是hard link，又称为硬链接；另一种是symbolic link，又称为符号链接。

-s：对源文件建立符号链接，而非硬链接
-f：强制创建链接，即使目标文件已经存在
```





