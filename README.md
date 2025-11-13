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
### Step1
GEO 检索: ((Alzheimer's disease AND Hippocampus) AND "Homo sapiens"[porgn:__txid9606])

下载文件 (summary)

交付 LLM 初步筛选合适的样本 (筛选条件:source: AD患者死后海马体 + date type: 转录组)

纳入数据集:
- GSE173955
- GSE184942
- GSE263319
- GSE84422
- GSE48350
- GSE5281
