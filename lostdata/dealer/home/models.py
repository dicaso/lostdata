#!/usr/bin/env python
from bidali import LSD
import pandas as pd, numpy as np

@LSD.storeDatasetLocally
def get_THMYCN():
    """
    Source: ~/LSData/private/2015_Hyperplasia_ABC/uniData.tx
    Source: ~/LSData/private/2015_Hyperplasia_ABC/SampleAnnotation.txt

    TODO refactor using retro DEA function
    """
    from bidali.retro import importr, ro, pandas2ri, base
    #TH-MYCN mouse model incorporation
    limma = importr('limma')
    prepcore = importr('preprocessCore')
    thannot = pd.read_table(LSD.datadir+'2015_Hyperplasia_ABC/SampleAnnotation.txt')
    thannot.group = thannot.group.apply(lambda x: x.replace('.','_'))
    thdata = pd.read_table(LSD.datadir+'2015_Hyperplasia_ABC/uniData.txt',sep=' ')
    thdata.columns = thannot.group + (['_r'+str(i+1) for i in range(4)]*3*2)
    genotype = thannot.genotype
    age = thannot.age
    genotype_r = ro.r.relevel(ro.FactorVector(genotype),ref='WT')
    thdatanorm_r = prepcore.normalize_quantiles(base.as_matrix(thdata))
    thdatanorm_r.colnames = ro.StrVector(thdata.columns)
    thdatanorm_r.rownames = ro.StrVector(thdata.index)
    metadata_r = base.data_frame(genotype=genotype_r,age=age)
    designSamples_r = ro.r['model.matrix'](ro.r.formula('~genotype:age'),data=metadata_r)
    fit_r = limma.lmFit(thdatanorm_r,designSamples_r)
    fit_r = limma.eBayes(fit_r)
    coefficients_r = fit_r.rx2('coefficients')
    thcoeffs = pd.DataFrame(
        pandas2ri.ri2py(coefficients_r),
        index=coefficients_r.rownames,
        columns=coefficients_r.colnames
    )
    thcoeffs['lineardiff'] = thcoeffs['genotypeTG:age']-thcoeffs['genotypeWT:age']

    return LSD.Dataset(exprData=thdata,
                       metadata=thannot,
                       lmCoeffs=thcoeffs)
