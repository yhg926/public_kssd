# Accuracy of resemblance estimation
---
###### By Huiguang Yi
Email: yhg926@gmail.com
###### 2019-08-08
---
## 0. Introduction
This is the workflow used to generate the results in Figure 1. Accuracy of resemblance estimation; 

Softwares versions in this experiment:
kssd version 1.0, bindash version 0.2.1 and Mash version 2.0.   

This experiment performed under linux system, 32G and 12-cores machine.   

## 1. Dataset
### 1.1 Testing datasets
The folder `dist1_30` includes AE016877.fasta and its 300 mutants with mutation rates range from 0.001 ~ 0.300; The folder `dist31_60` includes AE016877.fasta and its 300 mutants with mutation rates range from 0.301 ~ 0.600. The mutation rate of each mutant showed as the decimal between 'AE016877.' and '.fasta', for example: AE016877.0.128.fasta has a mutation rate of 0.128 from AE016877.fasta. These mutants were generated using the script `fasta_mut.pl`.
### 1.2 kssd_shuf_files
The .shuf files used by kssd for sketching.
### 1.3 Sketch-sizes file
bindash.sizes.txt and mash.size.txt listed the sketch-sizes for bindash and mash, respectively. 
## 2. Methods  
### 2.1 kssd
#### 2.1.1 sketching 

We first sketched `dist1_30` and `dist31_60` by dimensionality-reduction levels `$l$` = {4, 3, 2, 1, 0} and `$k$`= {8,10}, where `$k$` is half length of k-mer length in kssd (k=16,20). The folder `kssd_shuf_files` contains all the '.shuf' files used in the kssd resemblance estimation (manuscript Figure 1) except `$l$`= 0, which do not need a '.shuf' file (no dimensionality-reduction); File names follow this format: L`$l$`K`$k$`_`$group$`.shuf, where `$l$`, `$k$` and `$group$` are the arguments of dimensionality-reduction level, k-mer length and mutants groups, respectively. Here is an example for sketching: 
```
kssd dist -L kssd_shuf_files/L3K8_dist1_30.shuf -r dist1_30 -o L3K8_dist1_30
```
You can sketching using other .shuf files to get sketches with variated parameters.

For `$k$`=10 && `$l$` < 2, the k-mer substring space overflows this kssd implementation; then kssd need to be re-compiled; go to the kssd installation directory and type:
```
make strange
```
and run using the 'strange' kssd implementation, for example:
```
kssd dist -L 0 -k 10 -r dist1_30 -o L0K10_dist1_30
```

but do not use this 'strange' kssd implementation unless necessary, it is inefficient.  

#### 2.1.2  Distance estimation
The pair-wise distance file `L3K8_dist1_30/distance` were generated as follow (use `L3K8_dist1_30` as an example here):
```
kssd dist -r L3K8_dist1_30/ref -o L3K8_dist1_30 L3K8_dist1_30/qry
```
distance file for all other sketches were generated similarly.

#### 2.1.3  Correlation with ground truth mutation rates
Use  `L3K8_dist1_30` as an example:
```
cd L3K8_dist1_30 &&
awk '$1 ~ /dist1_30\/AE016877.fasta/ && $2 !~ /dist1_30\/AE016877.fasta/' ./distance.out|cut -f2,5|grep -v 'inf'|sed 's/\(dist1_30\/AE016877\.\|\.fasta\)//g'|sort -k1 -g|Rscript -e ' v<-read.table(file("stdin"),sep="\t");cor(v[,1],v[,2])'

[1] 0.9958683
```
all other correlation-coefficients were obtained similarly.


### 2.2 bindash
#### 2.2.1 Sketching

```
for i in `cat bindash.size.txt`; do bindash sketch --kmerlen=16 --sketchsize64="$i" --outfname=bindash.dist1_30k16s"$i".sketch --listfname=bindash_dist1_30.list --nthreads=12; done
```
`bindash.size.txt` includes all sketch-sizes set for the option --sketchsize64; `bindash_dist1_30.list` contain fullpaths for sequences in `dist1_30`. This is an example sketching using bindash with k=16 on `dist1_30`; You can perform sketching similarly using `$k$`= {20,21}, and on folder `dist31_60`.

#### 2.2.1 Distance estimation
```
mkdir bindash_out &&
mv *.sketch bindash_out/ &&
cd bindash_out/ &&
for i in `find ./*.sketch`; do bindash dist $i $i >$i.distance;done
```
#### 2.2.2 Correlation with ground truth mutation rates
```
 for i in `find ./bindash.dist1_30k20s*.sketch.distance`; do j=`awk '$1~/\/AE016877.fasta/ && $2!~/\/AE016877.fasta/' "$i" |cut -f2,3|sed 's/\(\.\.\/dist1_30\/AE016877\.\|.fasta\)//g' | sort -k1 -g |Rscript -e ' v<-read.table(file("stdin"),sep="\t");cor(v[,1],v[,2])'`;echo $i" "$j; done|sed 's/\(\.\/bindash\.dist1_30k20s\|\.sketch\.distance \[1\]\)//g' |sort -k1 -g
```
This gave Correlation-coefficients of bindash distance (of all sketch size) with ground truth, using `$k$` = 20, on folder `dist1_30`. You can do similarly using `$k$`={16,21}, and on folder `dist31_60`.  

### 2.3 Mash
#### 2.3.1 Sketching
```
for j in `cat mash.size.txt|grep -v '5236120'`; do mash sketch -k 20 -s $j ../dist31_60/*.fasta -o mash_k20s"$j"d31_60 -p 12;done
```
This is an example for mash sketching using `dist31_60`, k=20 and all sketch-sizes in mash.size.txt except '5236120'; sketch-size '5236120' caused running error in my machine, but you can try it on your machine. And you can perform sketching similarly using `$k$`= {16,21}, and on folder `dist1_30`.

#### 2.3.2 Distance estimation
```
mkdir mash_out &&
mv *.msh mash_out/ &&
cd mash_out &&
for i in `find ./*.msh`; do mash dist $i $i > $i.distance;done
```


#### 2.3.3 Correlation with ground truth mutation rates
```
 for i in `find *.distance|grep -v '5236120'`; do j=`awk '$1 ~ /\/AE016877.fasta/ && $2 !~ /\/AE016877.fasta/' $i | cut -f2,3 | sed 's/\(\.\.\/dist.*AE016877\.\|\.fasta\)//g' |sort -k1 -g|Rscript -e ' v<-read.table(file("stdin"),sep="\t");cor(v[,1],v[,2])' `;echo $i" "$j;done |grep 'd31_60'|grep 'mash_k21s'|sed 's/\(mash_k21s\|d31_60\.msh\.distance \[1\]\)//g' |sort -k1 -g
```
This gave Correlation-coefficients of mash distance (of all sketch-sizes except 5236120) with ground truth, using `$k$` = 21, on folder `dist31_60`.  You can do similarly using `$k$`={16,20}, and on folder `dist1_30`.
