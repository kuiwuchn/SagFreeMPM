Sagfree elastic MPM Demo
================
![Screenshot](demo.gif)
This is a Python implementation of initializating sag-free MPM simulation (https://graphics.cs.utah.edu/research/projects/sag-free-simulations/).

Dependencies
--------------------
This demo depends on [Taichi](https://github.com/taichi-dev/taichi). You may follow its [documentation](https://docs.taichi.graphics/) for installation, or simply type `pip install taichi` under a Python 3 environment.
This demo is tested under Taichi 1.0.4.

Run the Demo
--------------------
`python SagfreeElasticBeam.py`


BibTex Citation
----------------------
```
@article{Hsu2022,
   author       = {Jerry Hsu and Nghia Truong and Cem Yuksel and Kui Wu},
   title        = {A General Two-Stage Initialization for Sag-Free Deformable Simulations},
   journal      = {ACM Transactions on Graphics (Proceedings of SIGGRAPH 2022)},
   year         = {2022},
   month        = {07},
   volume       = {41},
   number       = {4},
   pages        = {64:1--64:13},
   articleno    = {64},
   numpages     = {13},
   location     = {Vancouver, Canada},
   url          = {https://doi.org/10.1145/3528223.3530165},
   doi          = {10.1145/3528223.3530165},
   issn         = {0730-0301},
   publisher    = {ACM Press},
   address      = {New York, NY, USA},
}
```