#1.Define clusters
less gencode.v42.primary_assembly.annotation.gtf |awk '$3=="exon"'|cut -f1,4,5,7,9 |cut -f1,2 -d";"|awk '{print$1,$2-1,$3,$8,$6,$4}'|sed -e 's/ /\t/g' | sed -e 's/[\";]//g' |sortBed >hg38exon.bed
ls -d *[0-9]|while read id ;do less ${id}/editing_sorted.txt |awk -F"\t" '{print$1,$2,$4,$9,$15,$16,$8}'|cut -f1,2,3,4,5,6,7 -d" "|awk 'BEGIN{OFS="\t"} $7=="CT"'|awk 'BEGIN{OFS="\t"} {if($3=="1"){print$1,$2-1,$2,$7,$4,"+",$5,$6;}else {if($3=="0"){print$1,$2-1,$2,$7,$4,"-",$5,$6;} else {print$1,$2-1,$2,$7,$4,".",$5,$6;}}}' OFS="\t"|sortBed|bedtools intersect -a - -b hg38exon.bed -wb -s|awk 'BEGIN{OFS="\t"} {print$12,$1,$3,$5,$6,$7,$8,$13}' |GenoTransPosShuttle.pl -o geno2trans -g gencode.v42.primary_assembly.annotation.gpe -i -|awk '{print$1,$9-1,$9,$2,$3,$5,$4,$6,$7,$8}' |sort -k1,1 -k2,2n -k3,3n |bedtools cluster -i - -d 100 -s|bedtools groupby -g 11 -c 1,1,10,4,5,5,6,7,8,9 -o count,mode,distinct,mode,min,max,mode,collapse,mode,mode|awk 'BEGIN{OFS="\t"} $2>9{print}'|sort -k4,4|bedtools groupby -g 4 -c 2,3,5,6,7,8,10,11 -o max,distinct,mode,min,max,mode,distinct,distinct|awk 'BEGIN{OFS="\t"} {print$4,$5-1,$6,$1,$2,$7,$8,$9}'|sortBed> newcluster/${id##*/}_cluster10.bed;done

#2.Transfer into transcriptome coord, choose the shortest length with the most mutations
bedtools intersect -a exp_editing_sorted.bed -b specific_cluster.bed -s >specific_edit_in_clu.bed
less hg38exon.bed |grep -E '^chr[0-9]*|[X-Y]\t'|bedtools intersect -a specific_edit_in_clu.bed -b - -s -wb | sort -k10|bedtools groupby -g 10 -c 1,1,2,3,6,11 -o mode,count,min,max,mode,distinct|awk 'BEGIN{OFS="\t"} {print$1,$2,$4+1,$3,$6,$7"\n"$1,$2,$5,$3,$6,$7}'|GenoTransPosShuttle.pl -o geno2trans -g gencode.v42.primary_assembly.annotation.gpe -i - |bedtools groupby -g 1 -c 2,3,3,4,5,6,7,7 -o mode,min,max,mode,mode,mode,min,max|awk 'BEGIN{OFS="\t"} {print$1,$2,$3-1,$4,$5,$6,$7,$8-1,$9,$9-$8+1}'|sort -k7 >transcript.txt
less transcript.txt |awk 'BEGIN{OFS="\t";max[key]=0} {key=$7;} {if (max[key]<$5||(max[key]==$5&&min[key]>$10)) {max[key]=$5;min[key]=$10;line[key]=$0;}} END{for(key in line) {print line[key]}}'|awk 'BEGIN{OFS="\t"} NR==FNR{a[$4];next} {if($7 in a) {print$1,$8,$9,$10,$2,$3,$4,$5,$6}}' specific_cluster.bed - >specific_trans_clu.bed

#Extract transcript reads:
gffread gencode.v42.primary_assembly.annotation.gtf -g GRCh38.primary_assembly.genome.fa -w gencode.v42.transcript.fa
less specific_trans_clu.bed|awk 'BEGIN{OFS="\t"} {print$1,$2,$3,$5"|"$6"|"$7"|"$8"|"$9"|"$4}'|bedtools getfasta -fi gencode.v42.transcript.fa -bed - -name -fo specific_trans_clu.fa

#3.Choose background
perl filter_main_chr.pl gencode.v42.3UTR.bed  main_background_3utr.bed
less main_background_3utr.bed |sortBed > main_background_3utr_sort.bed

bedtools intersect -a main_background_3utr.bed -b polyA_all_editing_sorted.bed -s -v|bedtools intersect -a - -b reditool_rmsk_Alu.bed -s -wa|sort -k5|awk '{key = $5;value = $7;if (value < min[key] || min[key] == "") {min[key] = value;line[key] = $0;}}END{for (key in line) {print line[key];}}'|awk '$7>10&&$7<8000'|shuf -n 685 >back_alu.bed
less back_alu.bed |cut -f5|sort -u >a.txt
bedtools intersect -a main_background_3utr.bed -b polyA_all_editing_sorted.bed -s -v|bedtools intersect -a - -b reditool_rmsk_Alu.bed -s -v |sort -k5|awk '{key = $5;value = $7;if (value < min[key] || min[key] == "") {min[key] = value;line[key] = $0;}}END{for (key in line) {print line[key];}}'|grep -v a.txt|awk '$7>10&&$7<8000'|shuf -n 845 >back_nonalu.bed
cat back_alu.bed back_nonalu.bed |sortBed >back.bed

bedtools intersect -a main_background_3utr.bed -b reditool_rmsk_Alu.bed -s -wa|sort -k5|bedtools intersect -a - -b polyA_all_editing_sorted.bed -s -v|sort -k5|awk '{key = $5;value = $7;if (value < min[key] || min[key] == "") {min[key] = value;line[key] = $0;}}END{for (key in line) {print line[key];}}'|sort -k7,7n >min_alu_noedit_sort.bed

less alu.fa|grep '>'|sed 's/[>:]/\t/g'|cut -f2,4|sort -k1,1n >cluster_alu_length_sort.bed

bedtools intersect -a reditool_rmsk_Alu.bed -b min_alu_noedit_sort.bed -s -wb|bedtools groupby -g 11 -c 7,8,9,10,12,13,2,3 -o mode,mode,mode,mode,mode,mode,first,first|sort -k7,7n|awk 'BEGIN{OFS="\t";i=1}  NR==FNR{a[NR]=$1;b[NR]=$2;next} {if(i<=685 && $7>=a[i]) {{print $0,i,a[i],b[i]};i+=1}}' cluster_alu_length_sort.bed - |awk 'BEGIN{OFS="\t"} {l=$11;sta=$8;end=$8+l} {if(end>$4) {print$2,$4-l,$4,$5,$1,$6,$7,$12,$11;} else{print$2,sta,sta+l,$5,$1,$6,$7,$12,$11;}}'>back_alu_length.bed
#control  the length of the selections and clusters' length are exactly the same :ensure that the length of the utr is longer than the length of the cluster, subtract the interval length from the start site of Alu, if within the interval, take it, if outside the interval, use the end of the interval to intercept the interval

less nonalu.fa|grep '>'|sed 's/[>:]/\t/g'|cut -f2,4|sort -k1,1n >cluster_nonalu_length_sort.bed
bedtools intersect -a main_background_3utr.bed -b polyA_all_editing_sorted.bed -s -v|bedtools intersect -a - -b back_alu_length.bed -s -v|sort -k5|awk '{key = $5;value = $7;if (value < min[key] || min[key] == "") {min[key] = value;line[key] = $0;}}END{for (key in line) {print line[key];}}'|bedtools intersect -a - -b reditool_rmsk_Alu.bed -s -v|sort -k7,7n > min_nonalu_noedi_sort.bed
awk 'BEGIN{OFS="\t";i=1}  NR==FNR{a[NR]=$1;b[NR]=$2;next} {if(i<=845 && $7>=a[i]) {{print $0,i,b[i],a[i]};i+=1}}' cluster_nonalu_length_sort.bed  min_nonalu_noedi_sort.bed|awk 'BEGIN{OFS="\t"} {print $1,$2,$2+$10,$4,$5,$6,$7,$9,$10}' >back_nonalu_length.bed

ls back*length.bed|while read id;do bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed ${id} -name -s -fo ${id%.*}.fa;done
cat back*length.fa >back_length.fa
cat back*length.bed|sortBed > back_length_sort.bed

#4.Secondary structure prediction
#RNAstructure
#exp:
less alu.fa |awk 'BEGIN{c=0;} {a[c++]=$0;} c==2{{system("echo \""a[0]"\n"a[1]"\" | Fold-smp -mfe - tmp.ct"); system("cat tmp.ct >>rs_exp_alu.ct")}; c=0;}'
rm tmp.ct
less nonalu.fa |awk 'BEGIN{c=0;} {a[c++]=$0;} c==2{{system("echo \""a[0]"\n"a[1]"\" | Fold-smp -mfe - tmp.ct"); system("cat tmp.ct >>rs_exp_nonalu.ct")}; c=0;}'
rm tmp.ct

#back:
less back_alu_length.fa |awk 'BEGIN{c=0;} {a[c++]=$0;} c==2{{system("echo \""a[0]"\n"a[1]"\" | Fold-smp -mfe - tmp.ct"); system("cat tmp.ct >>rs_back_alu.ct")}; c=0;}'
rm tmp.ct
less back_nonalu_length.fa |awk 'BEGIN{c=0;} {a[c++]=$0;} c==2{{system("echo \""a[0]"\n"a[1]"\" | Fold-smp -mfe - tmp.ct"); system("cat tmp.ct >>rs_back_nonalu.ct")}; c=0;}'
rm tmp.ct

#RNAfold
#back：
RNAfold <back_alu_length.fa --noPS >rf_back_alu.db; less rf_back_alu.db|cut -f1 -d" "|awk 'BEGIN{c=0;} {a[c++]=$0;} c==3{system("echo \""a[0]"\n"a[1]"\n"a[2]"\" | dot2ct - rftmp.ct"); system("cat rftmp.ct >>rf_back_alu.ct"); c=0;}';rm rftmp.ct
RNAfold <back_nonalu_length.fa --noPS >rf_back_nonalu.db; less rf_back_nonalu.db|cut -f1 -d" "|awk 'BEGIN{c=0;} {a[c++]=$0;} c==3{system("echo \""a[0]"\n"a[1]"\n"a[2]"\" | dot2ct - rftmp.ct"); system("cat rftmp.ct >>rf_back_nonalu.ct"); c=0;}';rm rftmp.ct
#exp：
RNAfold <specific_trans_clu.fa --noPS >exp.db; less exp.db|cut -f1 -d" "|awk 'BEGIN{c=0;} {a[c++]=$0;} c==3{system("echo \""a[0]"\n"a[1]"\n"a[2]"\" | dot2ct - rftmp.ct"); system("cat rftmp.ct >> specific_rf.ct"); c=0;}';rm rftmp.ct

#with SHAPE-seq
#The shape file was split into each gene, and the corresponding sequence was extracted according to the coordinates of the cluster

less shape.out.txt |awk 'BEGIN{OFS="\t"} {name=$1;l=$2+3;} {for (i=4;i<=l;i++) {print(i-3,$i,name)}}' >all.txt
sed -i 's/NULL/-999/g' all.txt
less all.txt |cut -f3|sort -u >shape_genelist.txt 
awk 'BEGIN{OFS="\t"} NR==FNR{a[$1];next} {if ($1 in a) print}' shape_genelist.txt specific_trans_clu.bed | awk 'BEGIN{OFS="\t"} {n=$1;a=$2;b=$3;} {cmd="less all.txt|grep "n"|head -n "b" |tail -n "b-a"|cut -f2|nl - >"n".txt"; system(cmd);}'

ls ENST*|while read id;do a=`less ${id}|wc -l`;b=`less ${id}|grep -w '\-999'|wc -l`;if [ $a == $b ];then echo ${id}>>nosignal.txt;fi;done
grep -f nosignal.txt
#exclude sequences without any signals

less shape_cluster_genelist.txt|while read id;do RNAfold < ${id}.fa --shape ${id}.txt  --noPS > ${id}.db;done
ls *.db|while read id; do less ${id}|cut -f1 -d" "|dot2ct - ct/${id%.*}.ct ;done
less nosignal.txt |while read line;do rm ${line%.*}.ct;done
cat ct/ENST00000* >>rf_shape_exp.ct

#shape_back:
mkdir db ct
less back_shape_signal.txt|while read id;do RNAfold < ${id}.fa --shape ${id}.txt  --noPS > db/${id}.db;done
cd db
ls *.db|while read id; do less ${id}|cut -f1 -d" "|dot2ct - ../ct/${id%.*}.ct ;done
cd ../
cat ct/*.ct > rf_shape_back.ct

#run RNAstructure with SHAPE-seq:
awk 'BEGIN{OFS="\t"} NR==FNR{a[$1];next} {if ($1 in a) print}' shape_genelist.txt specific_trans_clu.bed|cut -f1 > shape_cluster_genelist.txt
less shape_cluster_genelist.txt |while read line;do less specific_trans_clu.bed |grep ${line}|awk 'BEGIN{OFS="\t"} {print$1,$2,$3,$5"|"$6"|"$7"|"$8"|"$9"|"$4}'|bedtools getfasta -fi gencode.v42.transcript.fa -bed - -name -fo ${line}.fa;done
less shape_cluster_genelist.txt|while read id;do Fold-smp -mfe -sh ${id}.txt ${id}.fa  ${id}.ct;done
less nosignal.txt |while read line;do rm ${line%.*}.ct;done


#back_shape:
awk 'BEGIN{OFS="\t"} NR==FNR{a[$1];next} {if($4 in a) print}' shape_genelist.txt back_length_sort.bed >back_shape.bed


less back_shape.bed|cut -f4 >back_shape_transid.txt
less shape.out.txt |grep -f back_shape_transid.txt > back_shape.txt
less back_shape.txt |awk 'BEGIN{OFS="\t"} {name=$1;l=$2+3;} {for (i=4;i<=l;i++) {print(i-3,$i,name)}}'|sed 's/NULL/-999/g' >back_all.txt
less back_shape.bed |awk 'BEGIN{OFS="\t"} {print$4,$1,$2+1"\n"$4,$1,$3}'|GenoTransPosShuttle.pl -o geno2trans -g gencode.v42.primary_assembly.annotation.gpe -i -|bedtools groupby -g 1 -c 4,4 -o min,max|awk 'BEGIN{OFS="\t"} {print$1,$2-1,$3}'|sortBed >back_shape_transcoord.bed

less back_shape_transcoord.bed| awk 'BEGIN{OFS="\t"} {n=$1;a=$2;b=$3;} {cmd="less back_all.txt|grep "n"|head -n "b" |tail -n "b-a"|cut -f2|nl - >"n".txt"; system(cmd);}'
ls ENST*|while read id;do a=`less ${id}|wc -l`;b=`less ${id}|grep -w '\-999'|wc -l`;if [ $a == $b ];then echo ${id}>>nosignal.txt;fi;done

less nosignal.txt |sed 's/.txt//g'|cat - back_shape_transid.txt |sort |uniq -c|awk '$1<2{print$2}' > back_shape_signal.txt


cd back_shape
less back_shape_transid.txt |while read line;do less back_shape.bed |grep ${line}|bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed - -s -fo ${line}.fa;done

mkdir ct
less back_shape_signal.txt|while read id;do Fold-smp -mfe -sh ${id}.txt ${id}.fa  ct/${id}.ct;done
cat ct/ENST*.ct > rs_shape_back.ct
