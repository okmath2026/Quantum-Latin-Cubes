# QLS(5) 和 QLS(6) 小阶基数谱实验摘要

日期：2026-06-11

## 已确认的严格构造

### QLS(5)

由 5 阶 Latin square

```text
0 1 2 3 4
1 0 3 4 2
2 3 4 0 1
3 4 1 2 0
4 2 0 1 3
```

在左上角 intercalate

```text
0 1
1 0
```

上做一次二维旋转，得到基数 `c=7` 的 `QLS(5)`。

`intercalate_census.py` 已穷举全部 161280 个 5 阶 Latin squares，结果是：

- 最大 raw intercalate 数为 4；
- 最大可同时安全旋转的 compatible intercalate 数为 1；
- 因而纯 intercalate 旋转路线在 `n=5` 只能给 `c=5,7`。

### QLS(6)

对 cyclic Latin square `L(r,c)=r+c mod 6`，其 9 个 intercalate 把全部 36 个格子分成 9 个 2x2 块。选择不同子集旋转并做精确正交验证，得到

```text
{6, 8, 10, 12, 14, 16, 18} ⊆ Spec(6).
```

注意：旋转全部 9 个 intercalate 时，原来的 6 个计算基向量全部消失，所以基数是 18，而不是 `6+18=24`。

由 row-quantum Latin rectangle 乘积构造 `6=2·3`，还可得到

```text
{6, 10, 12, 18, 20, 24, 30, 36} ⊆ Spec(6).
```

合并两条严格路线：

```text
{6, 8, 10, 12, 14, 16, 18, 20, 24, 30, 36} ⊆ Spec(6).
```

其中 `c=7` 由 Zhang-Wang-Ji Lemma 3.1 排除。

## 数值证书

当前 `qls5_certificates.npz` 可读证书覆盖：

```text
QLS(5): c = 5, 7, 12, 21, 24, 25.
```

其中：

- `c=5,7` 可由显式 intercalate/经典构造严格化；
- `c=12` 已提取为完全精确的整数向量构造，见 `exact_qls5_c12_construction.py` 和 `qls5_c12_exact_note.md`；
- `c=21,24,25` 残差约 `2.2e-16`；

## 关键开放低阶目标

### QLS(5)

优先目标：

```text
c = 8, 9, 10, 11, 13, ..., 20, 22, 23.
```

尤其是 `c=8`。若不存在，说明“所有 n≥5 满谱除 n+1”的大猜想为假；若存在，则需要非 intercalate 型局部结构。

新增确定性筛选结果：`c8_pattern_search.py` 固定首行为计算基，枚举了 `c=8` 情形下 3 个新向量类的填格模式。运行结果：

```text
nodes = 3821945
leaves = 223104
unique_patterns = 34560
subspace_tests = 34560
feasible_count = 0
```

每个模式都会导出 3 个新向量的零坐标约束和相互正交图。所有 34560 个模式都落入以下三类线性代数障碍：

```text
same_2d_path:              8640
same_2d_triangle:         23040
two_2d_one_3d_triangle:    2880
```

解释：前两类都把三个新向量压进同一个二维坐标子空间，并强迫正交路径或正交三角；第三类强迫两个新向量成为同一二维子空间的一组正交基，第三个向量只能退化为计算基向量。

这强烈指向：

```text
8 ∉ Spec(5).
```

但要作为定理发表，还需要把三向量子空间筛选从随机线性代数检查改写成完全确定的有限情形证明。

详细证明笔记见 `qls5_c8_nonexistence_note.md`。

进一步整理：

- `c8_classification_certificate.py` 重新枚举并生成紧凑证书；
- `c8_classification_certificate.json` 给出分类计数与 SHA-256 摘要；
- `qls5_c8_theorem_draft.md` 已写成论文式定理证明草稿；
- `c8_verifier_audit.md` 记录验证器假设与无损剪枝。

证书摘要：

```text
unique_patterns = 34560
unobstructed_patterns = 0
sha256 = 8837bc4dd3c92ab424d30c4f1be7f20f235de8e87cee6563d6c202a609ed1a75
```

进一步压缩：`c8_active_graph_types.py` 把 34560 个模式按新向量所在格子的二部图同构类型压缩为 5 类：

```text
count   row degrees     column degrees     obstruction(s)
2880    (2,2,0,0)       (2,2,0,0,0)        same_2d_path
9600    (2,2,2,0)       (2,2,2,0,0)        same_2d_path / same_2d_triangle
14400   (3,2,2,0)       (3,2,2,0,0)        same_2d_triangle / two_2d_one_3d_triangle
5760    (3,3,2,0)       (3,3,2,0,0)        same_2d_triangle
1920    (3,3,3,0)       (3,3,3,0,0)        same_2d_triangle
```

这给手写分类证明提供了更清楚的路线：先分类 active-cell 二部图，再处理五个骨架上的标签填入。

更细的图层分析：

- `c8_active_graph_feasibility.py` 只使用 active-line、旧标签可填、新标签可染色三个必要条件，得到 35 个 active 图类型；
- `c8_active_graph_comparison.py` 比较后发现，完整标签与支持非退化约束淘汰其中 30 类，只剩上述 5 类；
- 5 个幸存类型具有同一个平衡特征：非零行度多重集等于非零列度多重集，即 `(2,2)`, `(2,2,2)`, `(3,2,2)`, `(3,3,2)`, `(3,3,3)`。

新的手写证明目标：证明“若不强迫某个新向量退化为计算基向量，则 active 图必须满足上述 balanced-degree 条件”。

新增 `c8_local_set_filter.py`：只使用局部集合条件

```text
A_r = row r 中出现的旧标签集合
B_c = {c} ∪ column c 中首行以下出现的旧标签集合
若新标签出现在 (r,c)，则 |A_r ∪ B_c| ≤ 3
```

并把同一新标签多次出现造成的零集累积也计入。结果正好是：

```text
necessary active graph types: 35
survive local set condition: 5
```

这说明 35 到 5 的压缩不需要向量数值搜索，可以纯粹看旧标签集合与支持非退化条件。

### QLS(5), c=12 的显式构造

从数值证书中提取出了一个精确构造。令 `e_i` 为计算基，设

```text
a = (0,-1,0, 2, 1)
b = (0, 1,0, 1,-1)
c = (0, 1,0, 0, 1)
d = (0,-1,0, 2,-5)
h = (0, 2,0, 1, 0)
f = (0, 1,0,-2, 0)
g = (0, 1,0, 0,-1)
```

则归一化下面的数组即可得到基数 12 的 `QLS(5)`：

```text
e0  e1  e2  e3  e4
a   e0  b   c   e2
d   e2  a   e0  h
h   e4  e0  e2  f
e2  e3  c   g   e0
```

因此现在有严格结果：

```text
{5,7,12} ⊆ Spec(5),   6,8 ∉ Spec(5).
```

`c=9` 也做了初步候选探测：`c9_pattern_probe.py` 固定首行为计算基，枚举 4 个新向量类的规范化标签模式，并用线性代数尝试赋向量。第一轮 120 秒结果：

```text
nodes = 1490916
leaves = 251570
unique_patterns = 9224
subspace_tests = 9224
found = False
```

这不是完整证明，但它提示 `c=9` 也不像简单小支持结构可达。结合已有 `c=12` 数值证书，值得重点研究是否

```text
Spec(5) 的低段为 {5,7,12,...}
```

或者至少存在连续缺口 `{8,9,10,11}`。

### QLS(6)

优先目标：

```text
c = 9, 11, 13, 15, 17, 19, 21, 22, 23, 25, ..., 35.
```

最关键的是第一个奇数目标 `c=9`。纯 intercalate 旋转只能给偶数增量，因此 `c=9` 必须来自更深的局部量子结构。

新增随机筛选：`c9_random_pattern_search.py` 固定首行为计算基，对 `c=9` 情形下的 3 个新向量类做随机约束填格，并用三向量子空间筛选验证。第一轮结果：

```text
attempts = 30
nodes = 151430
found = False
```

解释：这只是随机搜索，不能作为非存在性证据；但它说明 `QLS(6), c=9` 若存在，构型不会像 `n=5,c=7` 或 cyclic intercalate 那样直接。

## 下一步

1. 对 `QLS(5), c=8` 建立“三行以上修改”的结构化搜索。
2. 对 `QLS(6), c=9` 做有限向量池和数值 stratum search。
3. 将 `QLS(5), c=12` 证书代数化，尝试提取可发表的显式矩阵。
4. 把“两行刚性定理”写成正式引理：若除两行外全经典，则只能出现 `n+2t` 型基数。
