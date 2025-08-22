import os
import sys
from datetime import datetime
import contextlib
import subprocess
from urllib.request import urlretrieve
from urllib.request import HTTPError
from pathlib import Path
import glob
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder

class SabDabRetriever:
  def __init__(self):
    self.summary_file_path = None
    self.structures_dir = None
