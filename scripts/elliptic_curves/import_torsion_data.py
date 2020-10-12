# -*- coding: utf-8 -*-
r""" Import torsion growth data (as computed by Enrique Gonzalez).

Initial version (Warwick November 2017)

2018: updated for postgres BUT not yet checked to see if the data
types are what the postgres table wants.

Additional data fields for each elliptic curve over Q

 'tor_degs': u'jsonb'    list of degrees in which torsion grows
 'tor_fields': u'jsonb', list of field labels in which torsison grows

 'tor_gro': u'jsonb',    dictionary, keys are field labels or poly
                         strings, values are structure constants of the larger torsion

NB We have hard-coded the maximum degree of number field for which we
currently have data (currently 7) in lmfdb/elliptic_curves/web_ec.py
since tor_degs and the other extra columns are in the extr table.  If
data for higher degrees is uploaded that will need to be changed.
"""
from __future__ import print_function
import os
from sage.all import ZZ, PolynomialRing, QQ, NumberField, Set, copy

from lmfdb import db

print("setting curves")
curves = db.ec_curves
fields = db.nf_fields

Qx = PolynomialRing(QQ,'x')

def str_to_list(s):
    """
    Input: string representing list of ints, e.g. '', '2', '2,4', '5,4,3,2,1'
    """
    s = s.replace("[","").replace("]","")
    if s == '':
        return []
    else:
        return [int(c) for c in s.split(",")]

def poly_to_str(f):
    return str(f.coefficients(sparse=False)).replace(" ", "")[1:-1]


def find_field(pol, verbose=False):
    """
    pol is a string holding a list of coefficients, constant first, 1 last, e.g. '-2,0,1'

    Looks up this defining polynomial kn LMFDB and returns its label, or None
    """
    coeffs = str_to_list(pol)
    deg = len(coeffs)-1
    if deg==2:
        c, b, a = coeffs
        d = ZZ(b*b-4*a*c).squarefree_part()
        D = d if (d-1)%4==0 else 4*d
        absD = D.abs()
        s = 0 if d<0 else 2
        return '2.{}.{}.1'.format(s,absD)

    from lmfdb.number_fields.number_field import poly_to_field_label
    poly = Qx(coeffs)
    Flabel = poly_to_field_label(poly)
    if Flabel is None:
        print("********* field with polynomial {} is not in the database!".format(poly))
        K = NumberField(poly, 'a')
        poly = K.optimized_representation()[0].defining_polynomial()
        print("********* using optimised polynomial {}".format(poly))
        return poly_to_str(poly)
    else:
        if verbose:
            print("{} has label {}".format(pol,Flabel))
        return Flabel

def get_degree(label_or_coeffs):
    """Input: string, either a number field label or a list of
    coefficients of a defining polynomial, e.g. '-2,0,2'
    """
    if '.' in label_or_coeffs:
        return int(label_or_coeffs.split(".")[0])
    else:
        return label_or_coeffs.count(',')


def read_line(line, debug=0):
    r"""Parses one line from input file.  Returns a label and a dict with
    keys tor_degs, tor_fields, tor_gro.

    Sample line: 14a1 [3,6][1,1,1] [2,6][2,-1,1]

    Fields: label (single field, Cremona label)

            1 or more items of the form TF (with no space between)
            with T =[n] or [m,n] and F a list of integers of length
            d+1>=3 containing the coefficients of a monic polynomial
            of degree d defining a number field (constant coeff
            first).

    Note: in each file d is fixed and contained in the filename
    (e.g. growth2.000000-399999) but can be recovered from any line
    from the length of the coefficient lists.

    """
    if debug: print("Parsing input line {}".format(line))
    fields = line.split()
    label = fields[0]
    # print(fields[1:])
    # print([s[1:-1].split("][") for s in fields[1:]])
    tordata = [[F,T] for T,F in [s[1:-1].split("][") for s in fields[1:]]]
    tor_gro = dict(tordata)
    tor_fields = [F for F,T in tordata]
    tor_degs = sorted(list(Set([F.count(",") for F in tor_fields])))
    data = {'label': label,
            'tor_degs': tor_degs,
            'tor_fields': tor_fields,
            'tor_gro': tor_gro,
            }
    if debug: print("label {}, data {}".format(label,data))
    return label, data

all_degrees = [2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21]
all_ranges = ["0-9999"] + ["{}0000-{}9999".format(k,k) for k in range(1,40)]

HOME = os.getenv("HOME")
ECDATA_DIR = os.path.join(HOME, "ecdata")
GROWTH_DIR = os.path.join(ECDATA_DIR, "growth")

def read_torsion_growth_data(base_path=GROWTH_DIR, degrees=all_degrees, ranges=all_ranges):
    tor_data = {}
    count = 0

    def merge_data(d1, d2):
        # d1 and d2 are both dicts with keys tor_degs, tor_fields,
        # tor_gro return a merged dict, concatenating the tor_degs and
        # tor_fields lists and merging the tor_gro dicts
        #
        tor_gro = copy(d1['tor_gro'])
        tor_gro.update(d2['tor_gro'])
        return {'tor_degs': d1['tor_degs']+d2['tor_degs'],
                'tor_fields': d1['tor_fields']+d2['tor_fields'],
                'tor_gro': tor_gro,
        }

             
    for d in degrees:
        for r in ranges:
            filename = os.path.join(base_path, str(d), "growth{}.{}".format(d,r))   
            with open(filename) as h:
                print("opened %s" % filename)
                for line in h:
                    count += 1
                    if count%10000==0:
                        print("read %s lines" % count)
                    label, data = read_line(line)
                    if label in tor_data:
                        tor_data[label] = merge_data(tor_data[label], data)
                    else:
                        tor_data[label] = data

    print("finished reading {} lines from file".format(count))
    return tor_data

def write_torsion_growth_data(tordata, base_path='', filename="growth_table"):
    f = os.path.join(base_path, filename)
    with open(f,'w') as h:
        print("opened {}".format(f))

        # print header
        h.write("label|torsion_growth\n")
        h.write("text|jsonb\n\n")

        count = 0
        for lab, dat in tordata.items():
            tor_gro = str(dat['tor_gro']).replace("'",'"')
            h.write("|".join([lab, tor_gro]) + "\n")
            count +=1
            if count%1000==0:
                print("{} lines written so far".format(count))
    print("{} lines written to {}".format(count, f))
