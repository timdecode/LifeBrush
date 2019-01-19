# LifeBrush
A toolkit and environment for painting agent-based simulations, for mesoscale molecular illustrative simulations. This toolkit was orignally described in a conference paper, "LifeBrush: Painting interactive agent-based simulations" by Timothy Davison, Faramarz Samavati and Christian Jacob (the bibtex is below).

This toolkit has been improved significantly since the LifeBrush paper to design mesoscale molecular illustrative simulations. It is no longer generally applicable to any agent-based system through a mapping function (though, technically, that is still possible, the system just isn't designed around that anymore). The element framework and the agent framework are now the same.

Some things are broken as a result of major refactorings:
* exemplar design through painting operations (manual placement and configuration still works)
* synthesizing attribute vectors, they are replaced with a single attribute _type_

In the case of  synthesizing attribute vectors, I have limited things to just agent-types. Managing the mapping from types to attributes in the new agent-element framework was messy, and fragile. This is less general, but much easier to handle from the code side. I would like to explore a better alternative, that supports the attribute vectors, but doesn't break carefully configured agents.


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
