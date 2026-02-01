# PCBG
# One-Step Pairwise Constrained Multi-View Clustering in Linear Time

This repository is the official code implementation for the paper **"One-Step Pairwise Constrained Multi-View Clustering in Linear Time"**, corresponding to the work published in **IEEE Transactions on Knowledge and Data Engineering (TKDE)**.

The code implements our proposed **Pairwise Constrained Bipartite Graph (PCBG)** one-step multi-view clustering method, which achieves **linear time complexity** with respect to the number of samples.

---

## 📄 Paper Information

- **Title**: One-Step Pairwise Constrained Multi-View Clustering in Linear Time
- **Authors**: Wenjun Yu, Hong Tao, Chenping Hou
- **Journal**: IEEE Transactions on Knowledge and Data Engineering (TKDE)
- **Year**: 2025
- **Volume / Issue**: 37(9)
- **Pages**: 5523–5537
- **DOI**: [`10.1109/TKDE.2025.3579388`](https://doi.org/10.1109/TKDE.2025.3579388)
- **IEEE Xplore Link**: <https://ieeexplore.ieee.org/abstract/document/11033216>

---

## 🌟 Method Overview (PCBG Summary)

To address the issues in traditional multi-view clustering, such as:
- Requiring **multi-stage processing** (e.g., first constructing graphs then clustering),
- Difficulty in efficiently utilizing **pairwise constraint information** (must-link / cannot-link), and
- **Excessively high time complexity** for large-scale data,

This paper proposes the **Pairwise Constrained Bipartite Graph (PCBG)** learning method, achieving:

1. **One-Step Multi-View Clustering**: Within a unified optimization framework, simultaneously learning:
   - A consistent representation across multiple views
   - The propagation and encoding of pairwise constraints
   - A bipartite graph structure for clustering

2. **Linear Time Complexity**: The overall time complexity is linear with respect to the number of samples, making it more suitable for large-scale multi-view data scenarios.

3. **Pairwise Constraint Friendly**: Explicitly integrating must-link / cannot-link constraints during the construction of the bipartite graph, enabling efficient utilization of supervisory information.


## Method Framework Diagram

![PCBG Framework](figs/pcbg_framework.png)

*Figure: Schematic diagram of the proposed Pairwise Constrained Bipartite Graph (PCBG) one-step multi-view clustering framework.*


If you find our work useful for your research, please consider citing our paper:

```bibtex
@ARTICLE{11033216,
  author={Yu, Wenjun and Tao, Hong and Hou, Chenping},
  journal={IEEE Transactions on Knowledge and Data Engineering}, 
  title={One-Step Pairwise Constrained Multi-View Clustering in Linear Time}, 
  year={2025},
  volume={37},
  number={9},
  pages={5523-5537}
}
