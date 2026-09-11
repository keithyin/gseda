
分析比对到 reference 某个区域的 dw，ar 特征

分析流程：
1. 有两个 bam 输入 ，一个作为 对照组，一个作为实验组。
2. 将这两个 bam 比对到 reference ，使用 gsmm2 
3. 抽取比对到 reference 感兴趣的 region 上的 query 序列，query的base有其对应的 dw 和 ar 信息
4. 然后画出差异分布。 （dw：dwell time，是碱基信号的持续时间）（ar：arrival time，是从一个碱基转移到另一个碱基所持续的时间）
    对于 dw 来说，需要四个子图，分别对应 A、C、G、T，每个子图上有 实验组 和 对照组的 dw 分布图
    对于 ar 来说，16个子图，分别对应 A->A, A->C, A->G, A->T, .... 的 ar 
5. 由于两个bam的个数差异大，所以在画分布的时候，不要用 count，要用 frequency。画 hisgram 图。

注意：
1. 比对bam，是以 reference 方向为正方向，但是分析的过程中要以 query 的方向作为正
2. bam 中存在 ar，dw 字段，可以从中提取 ar，dw 信息，两个的长度均和 query base 的长度一致。在比对时，qeury 反向互补的时候，dw 和 ar 也被反向了
3. 关于 dw，以 query 为正向解释：每个 base 有一个对应的 dw-value，该dw-value 表示该 base 在测序时，信号持续了多长时间
4. 关于 ar，以 query 为正向解释：每个 base 有一个对应的 ar-value，该 ar-value 表示 前一个 base 信号结束 到 当前 base 信号开始，持续了多长时间