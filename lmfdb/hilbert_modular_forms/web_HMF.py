# -*- coding: utf-8 -*-
import re
import tempfile
import os
from pymongo import ASCENDING, DESCENDING
from flask import url_for, make_response
import lmfdb.base
from lmfdb.utils import comma, make_logger, web_latex, encode_plot

import sage.all
from sage.all import latex, matrix, ZZ, QQ

logger = make_logger("hmf")

hmf_forms = None
hmf_fields = None
nf_fields = None

def db_hmf_forms():
    global hmf_forms
    if hmf_forms is None:
        hmf_forms = lmfdb.base.getDBConnection().hmfs.forms
    return hmf_forms

def db_hmf_fields():
    global hmf_fields
    if hmf_fields is None:
        hmf_fields = lmfdb.base.getDBConnection().hmfs.fields
    return hmf_fields

def db_nf_fields():
    global nf_fields
    if nf_fields is None:
        nf_fields = lmfdb.base.getDBConnection().numberfields.fields
    return nf_fields

def construct_full_label(field_label, weight, level_label, label_suffix):
    if all([w==2 for w in weight]):           # Parellel weight 2
        weight_label = ''
    elif all([w==weight[0] for w in weight]): # Parellel weight
        weight_label = str(weight[0]) + '-'
    else:                                     # non-parallel weight
        weight_label = str(weight) + '-'
    return field_label + '-' + weight_label + level_label + '-' + label_suffix

class WebHMF(object):
    """
    Class for an Hilbert Modular Newform
    """
    def __init__(self, dbdata):
        """
        Arguments:

            - dbdata: the data from the database
        """
        logger.debug("Constructing an instance of WebHMF class")
        self.__dict__.update(dbdata)
        # All other fields are handled here
        self.make_form()

    @staticmethod
    def by_label(label):
        """
        Searches for a specific Hilbert newforms in the forms
        collection by its label.
        """
        data = db_ec().find_one({"label" : label})

        if data:
            return WebHMF(data)
        return "Hilbert newform %s not found" % label # caller must catch this and raise an error

    @staticmethod
    def from_data_string(label, L):
        """Takes an input line L from a raw data file and constructs the
        associated HMF object with given base field.

        String sample:
        <[31, 31, w + 12], "a", [-3, -2, 2, 4, -4, ...]>,
        """
        if isinstance(label, str):
            data['field_label'] = label
            F = HilbertNumberField(label)
            if "not found" in F:
                raise ValueError("No Hilbert nmber field with label %s is in the database" % label)
        elif label == None:
            raise ValueError("Must specify a valid field label")
        else: # we were passed a HilbertNumberField already
            F = label
            data['field_label'] = F.label

        # The level

        i = L.find('[')
        j = L.find(']')
        data['level_ideal'] = L[i,j+1].replace(" ","")
        N, n, alpha = data['level_ideal'][1:-1].split(',')
        data['level_norm'] = int(N)
        level = str2ideal(F,data['level_ideal'])
        data['level_label'] = F.ideal_label(level)

        # The weight

        data['parallel_weight'] = int(2)
        data['weight'] = str([data['parallel_weight']] * F.degree())

        # The label

        i = L.find('"')
        j = L.find('"', i+1)
        data['label_suffix'] = L[i,j+1].replace(" ","")
        data['label'] = construct_full_label(field_label,
                                             data['weight'],
                                             data['level_label'],
                                             data['label_suffix'])
        data['short_label'] = level_label + '-' + label_suffix

        # The hecke polynomial and degree

        if 'x' in L:
            # non-rational
            i = L.find("x")
            j = L.find(i+1,",")
            data['hecke_polynomial'] = L[i:j]
            data['dimension'] = int(1)
            x = polygen(QQ)
            hpol = x.parent()(str(pol))
            data['dimension'] = int(hpol.degree())
        else:
            # rational
            data['hecke_polynomial'] = 'x'
            data['dimension'] = int(1)

        i = L.rfind("[")
        j = L.rfind("]")
        data['hecke_eigenvalues'] = L[i+1:j].split(",")

        # Find (some of the) AL-eigenvalues

        BP = level.prime_factors()
        BP_indices = [F.ideal_index(P) for P in BP]
        BP_exponents = [level.valuation(P) fpr P in BP]
        AL_eigs = [int(data['hecke_eigenvalues'][i]) for i in BP_indices]
        if not all([(e==1 and eig in [-1,1]) or (eig==0)
                    for e,eig in zip(BP_exponents,AL_eigs)]):
            raise ValueError("Some bad AL-eigenvalues found")
        # NB the following will put 0 for the eigenvalue for primes
        # whose quare divides the level; this will need fixing later.
        data['AL_eigenvalues'] = [[F.primes[i],data['hecke_eigenvalues'][i]] for i in BP_indices]

        data['is_CM'] = '?'
        data['is_base_change'] = '?'
        data['AL_eigenvalues_fixed'] = None


    def save_to_db(self):
        pass

    def make_form(self):
        # To start with the data fields of self are just those from
        # the database.  We need to reformat these and compute some
        # further (easy) data about it.
        #
        pass
