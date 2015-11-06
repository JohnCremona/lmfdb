# -*- coding: utf-8 -*-
r""" Check presence of conjugates in table of Hilbert modular forms, and adds them if not present.

Assumes that the set of primes is a subset of the set of ideals.

Warning ; will not work with non-parallel weights because there is no
description of which weight corresponds to which embedding. (ordered by image of
generator in R ?)
Note that labels contain no information about weight.

Initial version (University of Warwick 2015) Aurel Page

"""

import sys
sys.path.append("../..");
import pymongo
from lmfdb import base
from lmfdb.website import dbport
from lmfdb.WebNumberField import WebNumberField
from lmfdb.hilbert_modular_forms.hilbert_field import (findvar, niceideals,
 conjideals, str2ideal, HilbertNumberField)
from pymongo import MongoClient

print "calling base._init()"
dbport=37010
base._init(dbport)
print "getting connection"
C = base.getDBConnection()
print "authenticating for the hmfs database"
C['hmfs'].authenticate('editor', '282a29103a17fbad')
print "authenticating for the numberfields database"
#FIXME: use readonly login/password
C['numberfields'].authenticate('editor', '282a29103a17fbad')
print "setting hmfs, fields and forms"
hmfs = C.hmfs
fields = hmfs.fields
forms = hmfs.forms
nfcurves = C.elliptic_curves.nfcurves

# Cache of WebNumberField and FieldData objects to avoid re-creation
WNFs = {}
Fdata = {}

def get_WNF(label, gen_name):
    if not label in WNFs:
        WNFs[label] = WebNumberField(label, gen_name=gen_name)
    return WNFs[label]

def get_Fdata(label):
    if not label in Fdata:
        Fdata[label] = fields.find_one({'label':label})
    return Fdata[label]

def checkprimes(label):
    Fdata = get_Fdata(label)
    gen_name = findvar(Fdata['ideals'])
    WebF = get_WNF(label, gen_name)
    F = WebF.K()
    ideals = niceideals(F, Fdata['ideals'])
    primes = niceideals(F, Fdata['primes'])
    F = HilbertNumberField(label)
    L = []
    for prhnf,prideal,prlabel in primes:
        ideal = F.ideal(prlabel)
        if ideal != prideal:
            L.append(prlabel)
    return L

def nautos(label):
    return len(get_WNF(label, 'a').K().automorphisms())

def fldlabel2conjdata(label):
    data = {}
    Fdata = get_Fdata(label)
    gen_name = findvar(Fdata['ideals'])
    WebF = get_WNF(label, gen_name)
    F = WebF.K()
    data['F'] = F
    auts = F.automorphisms()
    if len(auts) == 1: #no nontrivial automorphism, nothing to do
        return None
    auts = [g for g in auts if not g.is_identity()]
    data['auts'] = auts
    ideals = niceideals(F, Fdata['ideals'])
    data['ideals'] = ideals
    cideals = conjideals(ideals, auts)
    data['conjideals'] = cideals
    primes = niceideals(F, Fdata['primes'])
    data['primes'] = primes
    primeslabels = [prm[2] for prm in primes]
    cprimes = conjideals(primes, auts)
    cprimes = [[primeslabels.index(cprimes[(prm[2],ig)]) for prm in primes] for ig in range(len(auts))]
    data['conjprimes'] = cprimes
    return data

def conjstringideal(F,stridl,g):
    """Given a string representing an ideal of F and an automorpism g of F,
    return the string representing the conjugate ideal.
    """
    N,n,_,gen = str2ideal(F,stridl)
    return '[' + str(N) + ',' + str(n) + ',' + str(g(gen)) + ']'

def conjform_label(f, ig, cideals):
    level_label = cideals[(f['level_label'],ig)]
    short_label = level_label + '-' + f['label_suffix']
    return f['field_label'] + '-' + short_label

def conjform(f, g, ig, cideals, cprimes, F): #ig index of g in auts
    if f['is_base_change'][0:3] == 'yes':
        print("This form is a base-change.")
        #return None
        #if the form is a base-change, but not from Q,
        #we should still add its conjugate
    fg = copy(f)

    fg['level_label'] = cideals[(f['level_label'],ig)]
    fg['short_label'] = fg['level_label'] + '-' + fg['label_suffix']
    fg['label'] = fg['field_label'] + '-' + fg['short_label']

    fg['level_ideal'] = conjstringideal(F,f['level_ideal'],g)

    fg['AL_eigenvalues'] = [[conjstringideal(F,x[0],g),x[1]] for x in f['AL_eigenvalues']]

    H = f['hecke_eigenvalues']
    Hg = copy(f['hecke_eigenvalues'])
    fg['hecke_eigenvalues'] = Hg

    attained  = [False for i in range(len(H))]
    for i in range(len(H)):
        if cprimes[ig][i] < len(H):
            attained[cprimes[ig][i]] = True
    maxi = 0
    while maxi < len(H):
        if not attained[maxi]:
            break
        maxi += 1
    if maxi < len(H):
        print("truncating list of eigenvalues (missing conjugate prime)")
    del Hg[maxi:]

    for i in range(len(H)):
        if cprimes[ig][i] < maxi:
            Hg[cprimes[ig][i]] = H[i]

    del fg['_id']
    return fg

def checkadd_conj(label, min_level_norm=0, max_level_norm=None, fix=False, buildform=False):
    if fix:
        buildform = True
    count = 0
    countmiss = 0
    query = {}
    query['field_label'] = label
    query['level_norm'] = {'$gte' : int(min_level_norm)}
    if max_level_norm:
        query['level_norm']['$lte'] = int(max_level_norm)
    else:
        max_level_norm = oo
    ftoconj = forms.find(query)
    print("%s forms over %s to examine of level norm between %s and %s."
          % (ftoconj.count(),label,min_level_norm,max_level_norm))
    if ftoconj.count() == 0:
        return None
    print("Ideals precomputations...")
    data = fldlabel2conjdata(label)
    if data == None:
        print("No nontrival automorphisms!")
        return
    print("...done.\n")
    auts = data['auts']
    print("Applying %s non-trivial automorphisms..." % len(auts))
    cideals = data['conjideals']
    cprimes = data['conjprimes']
    F = data['F']
    for f in ftoconj:
        print("Testing form %s" % f['label'])
        for g in auts:
            ig = auts.index(g)
            fg_label = conjform_label(f, ig, cideals)
            fgdb = forms.find_one({'label':fg_label})
            if fgdb == None:
                print("conjugate not present")
                countmiss += 1
                if buildform:
                    fg = conjform(f, g, ig, cideals, cprimes, F)
                if fix:
                    if fg != None: #else: is a lift (self-conjugate), should have been detected
                        print("adding it : "+fg['label'])
                        forms.insert(fg)
                        count += 1
    print("\nMissing "+str(countmiss)+" conjugate forms (possibly counted multiple times if several nontrivial automorphisms).")
    print("Added "+str(count)+" new conjugate forms.")
    return None

def forms_equal(f,g):
    fH = f['hecke_eigenvalues']
    gH = g['hecke_eigenvalues']
    for i in range(min(len(fH),len(gH))):
        if fH[i] != gH[i]:
            return False
    return True

def check_multiplicity_one(label):
    F = HilbertNumberField(label)
    count = 0
    for N in F.ideals_iter():
        Lf = forms.find({'field_label':label, 'level_label':N['label']})
        Lf = [f for f in Lf]
        n = len(Lf)
        for i in range(n):
            for j in range(i+1,n):
                if forms_equal(Lf[i],Lf[j]):
                    count += 1
                    print "duplicates: "+Lf[i]['label']+" and "+Lf[j]['label']
    print("Found "+str(count)+" duplicate forms.")

def fix_data_fields(min_level_norm=0, max_level_norm=None, fix=False):
    r""" One-off utility to:
    1. add degree and disc fields for each Hilbert newform
    2. Change CM and base-change from "yes?" to "yes"
    """
    count = 0
    query = {}
    query['level_norm'] = {'$gte' : int(min_level_norm)}
    if max_level_norm:
        query['level_norm']['$lte'] = int(max_level_norm)
    else:
        max_level_norm = oo
    forms_to_fix = forms.find(query)
    print("%s forms to examine of level norm between %s and %s."
          % (forms_to_fix.count(),min_level_norm,max_level_norm))
    if forms_to_fix.count() == 0:
        return None
    for f in forms_to_fix:
        count = count+1
        if count%100==0: print("%s: %s" % (count, f['label']))
        fix_data = {}
        deg, r, disc, n = f['field_label'].split('.')
        fix_data['deg'] = int(deg)
        fix_data['disc'] = int(disc)
        if f['is_CM'] == 'yes?':
            fix_data['is_CM'] = 'yes'
        if f['is_base_change'] == 'yes?':
            fix_data['is_base_change'] = 'yes'
        #print("using fixed data %s for form %s" % (fix_data,f['label']))
        if fix:
            forms.update({'label': f['label']}, {"$set": fix_data}, upsert=True)

def fix_one_label(lab, reverse=False):
    r""" If lab has length 1 do nothing.  If it has length 2 increment the
    first letter (a to b to c to ... to z).  The lenths must be at
    most 2 and if =2 it must start with 'a'..'y' (these are all which
    were required).  If reverse==True the inverse operation is carried
    out (z to y to ... to c to b to a).
    """
    if len(lab)!=2:
        return lab
    else:
        if reverse:
            return chr(ord(lab[0])-int(1))+lab[1]
        else:
            return chr(ord(lab[0])+int(1))+lab[1]

def fix_labels(min_level_norm=0, max_level_norm=None, fix=False, reverse=False):
    r""" One-off utility to correct labels 'aa'->'ba'->'ca', ..., 'az'->'bz'->'cz'
    """
    count = 0
    query = {}
    query['level_norm'] = {'$gte' : int(min_level_norm)}
    if max_level_norm:
        query['level_norm']['$lte'] = int(max_level_norm)
    else:
        max_level_norm = oo
    forms_to_fix = forms.find(query)
    print("%s forms to examine of level norm between %s and %s."
          % (forms_to_fix.count(),min_level_norm,max_level_norm))
    if forms_to_fix.count() == 0:
        return None
    for f in forms_to_fix:
        count = count+1
        if count%100==0: print("%s: %s" % (count, f['label']))
        fix_data = {}
        lab = f['label_suffix']
        if len(lab)==1:
            continue
        if f['label'][-2:] != lab:
            print("Incorrect label_suffix %s in form %s" % (lab,f['label']))
            return
        oldlab = lab
        lab = fix_one_label(lab, reverse=reverse)
        fix_data['label_suffix'] = lab
        fix_data['label'] = f['label'].replace(oldlab,lab)
        fix_data['short_label'] = f['short_label'].replace(oldlab,lab)
        print("using fixed data %s for form %s" % (fix_data,f['label']))
        if fix:
            forms.update({'label': f['label']}, {"$set": fix_data}, upsert=True)

        # find associated elliptic curve and fix that too (where appropriate)
        if f['deg']==2 and f['dimension']==1:
            label = f['label']
            for e in nfcurves.find({'class_label':f['label']}):
                fix_data = {}
                fix_data['iso_label'] = lab
                fix_data['label'] = e['label'].replace(oldlab,lab)
                fix_data['short_label'] = e['short_label'].replace(oldlab,lab)
                fix_data['class_label'] = e['class_label'].replace(oldlab,lab)
                fix_data['short_class_label'] = e['short_class_label'].replace(oldlab,lab)
                print("using fixed data %s for curve %s" % (fix_data,e['label']))
                if fix:
                    nfcurves.update({'label': e['label']}, {"$set": fix_data}, upsert=True)
        else:
            print("No elliptic curve to fix")

def add_numeric_label_suffixes(min_level_norm=0, max_level_norm=None, fix=False):
    r""" One-off utility to add a numeric conversion of the letter-coded
    label suffixes 'a'->0', 'z'->25, 'ba'->26, etc. for sorting
    purposes.
    """
    from sage.databases.cremona import class_to_int
    count = 0
    query = {}
    query['level_norm'] = {'$gte' : int(min_level_norm)}
    if max_level_norm:
        query['level_norm']['$lte'] = int(max_level_norm)
    else:
        max_level_norm = oo
    forms_to_fix = forms.find(query)
    print("%s forms to examine of level norm between %s and %s."
          % (forms_to_fix.count(),min_level_norm,max_level_norm))
    for f in forms_to_fix:
        count = count+1
        if count%100==0: print("%s: %s" % (count, f['label']))
        fix_data = {}
        lab = f['label_suffix']
        fix_data['label_nsuffix'] = class_to_int(lab)
        #print("using fixed data %s for form %s" % (fix_data,f['label']))
        if fix:
            forms.update({'label': f['label']}, {"$set": fix_data}, upsert=True)


