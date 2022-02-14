# Evolution-of-Subfamily-I.1-Lipases-in-Pseudomonas-aeruginosa

此分析流程的目的在于重复一篇文章，以熟悉蛋白进化分析的一般流程
[《Evolution of Subfamily I.1 Lipases in Pseudomonas aeruginosa》](https://link.springer.com/article/10.1007/s00284-021-02589-4)

## 1. 数据下载

![](./IMG/Table1.png)

根据文章Table1提供的数据，从NCBI中选取相应的物种，下载其基因组

### 1.1 基因组数据（实验数据）

```bash
mkdir -p /mnt/d/project/Evolution/genome     #工作路径
```

### 1.2 脂肪酶序列数据（目标序列）

```bash
mkdir -p /mnt/d/project/Evolution/lipase
```

## 2. 检索不同假单胞菌中I.1脂肪酶的数量

### 2.1 BLAST工具下载

简介：

BLAST，全称Basic Local Alignment Search Tool，即"基于局部比对算法的搜索工具"，由Altschul等人于1990年发布。BLAST能够实现比较两段核酸或者蛋白序列之间的同源性的功能，它能够快速的找到两段序列之间的同源序列并对比对区域进行打分以确定同源性的高低。

BLAST的运行方式是先用目标序列建数据库（这种数据库称为database，里面的每一条序列称为subject），然后用待查的序列（称为query）在database中搜索，每一条query与database中的每一条subject都要进行双序列比对，从而得出全部比对结果。

BLAST是一个集成的程序包，可以实现五种可能的序列比对方式：

```
blastp：蛋白序列与蛋白库做比对，直接比对蛋白序列的同源性。

blastx：核酸序列对蛋白库的比对，先将核酸序列翻译成蛋白序列（根据相位可以翻译为6种可能的蛋白序列），然后再与蛋白库做比对。

blastn：核酸序列对核酸库的比对，直接比较核酸序列的同源性。

tblastn：蛋白序列对核酸库的比对，将库中的核酸翻译成蛋白序列，然后进行比对。

tblastx：核酸序列对核酸库在蛋白级别的比对，将库和待查序列都翻译成蛋白序列，然后对蛋白序列进行比对。
```
