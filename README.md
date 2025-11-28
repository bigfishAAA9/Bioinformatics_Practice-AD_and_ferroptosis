# Bioinformatics_Practice-AD_and_ferroptosis
## 了解背景知识
文献
- https://doi.org/10.1038/s41392-020-00428-9
- https://doi.org/10.1038/s41392-024-02071-0

机制图
<img width="2945" height="2216" alt="mermaid-20251112 203600" src="https://github.com/user-attachments/assets/f970133b-2e56-45df-8c81-96de734fdb5d" />

铁失衡热点脑区
| 脑区                                                     | 铁水平变化趋势 | 与病理关联                |
| ------------------------------------------------------ | ------- | -------------------- |
| **海马（hippocampus）**（尤其CA1）                             | ↑       | 铁沉积显著，与认知下降高度相关      |
| **颞叶皮层（temporal cortex / temporal lobe）**              | ↑       | 与Aβ沉积、tau病理、认知下降显著相关 |
| **额叶皮层（frontal cortex）**                               | ↑       | 与Aβ、tau病理正相关         |
| **新皮层（neocortex）**                                     | ↑       | 铁积累广泛                |
| **基底节（basal ganglia, 特别是putamen）**                     | ↑       | 部分研究发现高铁             |
| **苍白球（globus pallidus, GP）**                           | ↓       | 与PD不同趋势              |
| **扣带、顶叶、内嗅皮层（cingulate, parietal, entorhinal cortex）** | —       | 铁水平变化不显著             |

## 获取数据
GEO 检索: ((Alzheimer's disease AND Hippocampus) AND "Homo sapiens"[porgn:__txid9606])

人工筛选，得到准入的数据集表格，见 ./1step/收集表.xlsx

共有22个数据集，涵盖uid、GEOaccession、platform、sample_number(AD/NONAD)、sample_info五个维度信息

## 数据预处理
### Step1
对上述数据集再进行筛选，最后采用了10个数据集，表格可见'./1step/dataset_used.xlsx'，7个array数据集(共计243个样本)，3个RNA-Seq数据集(共计24个样本)

下载以上数据集，进行log(exp+1)、normalizeBetweenArrays处理

### Step2
使用sva包分别对array、RNA-Seq数据集进行去批次效应，以下展示array的结果，根据箱图，可以明显地看出，结果并不明显，在Step3DEG结果也可体现

#### 处理前
boxplot
  <img width="1587" height="1091" alt="image" src="https://github.com/user-attachments/assets/94d53ddf-f717-43b6-8447-601ce39bcc78" />
PCA
  批次分组
  <img width="1670" height="1155" alt="image" src="https://github.com/user-attachments/assets/5558eee4-f0e2-443d-a349-e919944f1ef4" />
  疾病状态分组
  <img width="1661" height="1154" alt="image" src="https://github.com/user-attachments/assets/ce94d24e-0f60-40ea-9e73-24c5a8829dea" />

#### 处理后
boxplot
  <img width="1597" height="1104" alt="image" src="https://github.com/user-attachments/assets/569d80e6-5562-4cda-b2ac-01a647301a73" />
PCA
  批次分组
  <img width="1722" height="1181" alt="image" src="https://github.com/user-attachments/assets/4c6f8fac-b224-4272-9cc5-f6d97b3749fe" />
  疾病状态分组
  <img width="1666" height="1153" alt="image" src="https://github.com/user-attachments/assets/32493177-d960-4a5f-9daa-555954c2151f" />


### Step3
使用limma包在array数据集中筛选DEGs，共得到705个差异基因( log2FC > 0.2 & padj < 0.05 )，与FRGs(264个基因)共有17个交集基因

volcano
  <img width="1670" height="1175" alt="image" src="https://github.com/user-attachments/assets/522121fc-2b90-4201-9c17-bd9c914cfbd8" />

veen
  <img width="715" height="517" alt="image" src="https://github.com/user-attachments/assets/25a58b97-0358-4817-a4b4-6decf63d4a03" />



