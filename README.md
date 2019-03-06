# LifeBrush
A toolkit and environment for painting agent-based simulations, for mesoscale molecular illustrative simulations. This toolkit was orignally described in a conference paper, "LifeBrush: Painting interactive agent-based simulations" by Timothy Davison, Faramarz Samavati and Christian Jacob (the bibtex is below).

This toolkit has been improved significantly since the LifeBrush paper to design mesoscale molecular illustrative simulations. It is no longer generally applicable to any agent-based system through a mapping function (though, technically, that is still possible, the system just isn't designed around that anymore). The element framework and the agent framework are now the same.

Some things are broken as a result of major refactorings and new ideas:
* exemplar design through painting operations (manual placement and configuration still works)
* synthesizing attribute vectors, they are replaced with a single attribute _type_

The synthesis of attribute vectors was fragile. It was easy to break agents during the mapping process. For now, synthesis copies the agent configuration directly from the exemplar to the output. It is still possible to reassign agent behaviors by creating a new exemplar with the required behavior, but this will reassign all of the agent attributes from that exemplar to the output agent. This method is much faster and cleaner, plus it isn't as fragile as the old idea. I'll probably keep it.

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
