#!/usr/bin/env python
import numpy as np
import yaml
import os
import ROOT

class Config:
  def __init__(self, path=None):
    self.doc = {}
    if path: self.load(path)

  def __str__(self):
    return yaml.dump(self.doc, indent=2, sort_keys = False)

  def load(self, path):
    if isinstance(path, str):
      if not os.path.exists(path):
        raise FileNotFoundError(f'Configuration file not found: {path}')

      with open(path) as f:
        docs = yaml.safe_load_all(f)
        for doc in docs:
          if isinstance(doc, dict):
            self._deep_merge_dict(self.doc, doc)
          else:
            print(f'Warning: skipping non-dict type doc in {path} YAML')

    elif isinstance(path, dict):
      self._deep_merge_dict(self.doc, path)

  def _deep_merge_dict(self, target, source):
    for key, value in source.items():
      if key in target and isinstance(target[key], dict) and isinstance(value, dict):
        self._deep_merge_dict(target[key], value)
      else:
        target[key] = value

  def get(self, key_path):
    keys = key_path.split('.')
    doc = self.doc
    for key in keys:
      if isinstance(doc, dict) and key in doc:
        doc = doc[key]
      else:
        return None
    return doc

class ConfigRENE(Config):
  def __init__(self, path):
    super().__init__(path)

  def getReactors(self, key=None):
    reactors = self.get('reactors')
    if key == None:
      return reactors
    else:
      retval = []
      for reactor in reactors:
        retval.append(reactor[key])
      return retval

  def getDetectors(self, key=None):
    detectors = self.get('detectors')
    if key == None:
      return detectors
    else:
      retval = []
      for detector in detectors:
        retval.append(detector[key])
      return retval

  def getBaselines(self, detName):
    ## Get the detector position x,y,z
    detNames = self.getDetectors('name')
    if detName not in detNames: return None
    detPoses = self.getDetectors('position')
    x, y, z = detPoses[detNames.index(detName)]

    baselines = []
    for srcX, srcY, srcZ in self.getReactors('position'):
      r = np.linalg.norm([x-srcX, y-srcY, z-srcZ])
      baselines.append(np.round(r, decimals=2))

    return baselines

def getFileAndObj(path):
  path = path.split(':', 1)
  fPath = path[0]
  hPath = path[1] if len(path) > 0 else None

  f = ROOT.TFile(fPath)
  h = f.Get(hPath)

  return f, h

if __name__ == '__main__':
    config = Config("config.yaml")
    print(config)
