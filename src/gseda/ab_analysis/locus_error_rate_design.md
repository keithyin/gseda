

实现一个 reference 不同位点的 准确率的 AB test 结果

输入：
    query：实验组数据，对照组数据。数据格式 bam、fastq 均可。（bam是 unmapped 的）
    ref：reference fasta
    rq-thr: read accuracy 的阈值。是一个list，可以填一个，可以填两个。如果是一个的话，实验组和对照组就共用一个。注意是 read-accuracy. 不是 phreq。 这个阈值会用在 gsmm2 比对时候

处理流程：
    query 都比对到 ref 上。使用 gsmm2
    然后生成的比对 bam 调用 gsetl ,gsetl 会输出一个 fact_aligned_bam_ref_locus_info.csv 文件，将这个文件处理一下，再增加一个每个位点准确率的数据。然后将实验组和对照组的两个表 join 起来 当成一个表。我需要这个新表做判断
    

脚本写在 @/root/projects/gsda/third_party/gseda/src/gseda/ab_analysis/locus_error_rate.py