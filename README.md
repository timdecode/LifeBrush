# LifeBrush
LifeBrush is a VR Unreal Engine 4 toolkit for creating interactive biomolecular illustrative simulations. It combines interactive painting of molecular patterns in VR, with an agent-based modeling framework for defining molecular interactions and behaviours. 

LifeBrush is fast. The generative brushes can synthesize thousands of molecules per second in volumes and on surfaces. The agent-based modeling framework is also fast, capable of simulating more than 10,000 agents at 90 fps in VR.

This toolkit was orignally described in a conference paper, "LifeBrush: Painting interactive agent-based simulations" by Timothy Davison, Faramarz Samavati and Christian Jacob (the bibtex is below). The version here on GitHub is an evolution of that paper, but it is a refactoring and drops the mapping idea. Now, the framework directly synthesizes agents into a simulation (this is simpler and cleaner, though less flexible). 

# License and Copyright 

All code, unless stated otherwise, is Copyright (c) 2019, Timothy Davison. All rights reserved.

All of my code is released under a MIT license. There are included source-codes released under their respective licenses.

**Please cite with:**
```
@inproceedings{davison2018lifebrush,
  title={LifeBrush: Painting interactive agent-based simulations},
  author={Davison, Timothy and Samavati, Faramarz and Jacob, Christian},
  booktitle={2018 International Conference on Cyberworlds (CW)},
  pages={17--24},
  year={2018},
  organization={IEEE}
}
```
