from ._processing import *
from ._fileIO import *
from ._sp_plot import * 
from ._utils import *
from ._enrichment import *
from ._cci import *
from ._dataset import Visium_sham, Xenium_sham, Xenium_iri, Xenium_iri_poly

def load_Visium_data():
    return Visium_sham.load_anndata()

def load_Xenium_data():
    return Xenium_sham.load_anndata()

def load_Xenium_iri_data():
    return Xenium_iri.load_anndata()

def load_Xenium_iri_poly():
    return Xenium_iri_poly.load_json()