
使用gsmm2 将 smc 输出的 bam/fastq 比对到 reference 上。

将 query 分成三组：
    1. query 在 reference 某个感兴趣的区间有错误（ins、del，mismatch）
    2. query 在 reference 某个感兴趣区间没有错误 （eq）
    3. query 没有比对到 reference 的某个感兴趣区间

基于 query 的名字，可以抽取到其对应的 channel，这样就会形成三组 channel。保留前两组，最后一组丢弃。
保留的两组 称为 错误组、正确组


对于 输入的 adapter.bam 进行处理，将其按照上述的 错误组、正确组的 channels 将该文件分成两个部分。通过 adapter.bam 中的 ch 字段进行关联
