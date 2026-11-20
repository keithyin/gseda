
要实现一个 STR、SSR、homopolymer 区域检测的功能，写成一个 python 脚本

1. 使用 trfrs --minscore 9 --str 进行粗筛
2. 粗筛结果进行处理，
    输出的结果不能有 overlap
    如果返回了类似 GCGC 的motif，那得归一化成 GC

3. 基于粗筛结果进行精筛，提供两个模式，一个 --str 模式，一个是 --ssr 模式
    --ssr 模式要求区域中必须是完全干净的 unit * repeats
    --str 直接用粗筛的结果即可

测试数据在 /data1/ccs_data/str-optimization/re-rerun/reference 目录下，这些都是 STR 相关的数据