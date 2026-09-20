import json
import numpy as np
import pandas as pd
from scipy import stats

from cross_validate_nasa import external_sets
from config import RESULTS_DIR

SITE = RESULTS_DIR

PANELS = {
    'GLDS-7 root defense-regulatory': {
        'AT2G35930':'PUB23','AT5G59820':'ZAT12/RHL41','AT1G19210':'ERF017',
        'AT1G13260':'RAV1','AT1G19180':'JAZ1','AT2G44840':'ERF13',
        'AT1G28370':'ERF11','AT5G47220':'ERF2','AT3G45640':'MPK3',
        'AT1G80840':'WRKY40','AT4G23810':'WRKY53','AT5G27420':'ATL31/CNI1'
    },
    'NASA cross-study core': {
        'AT1G74310':'HSP101','AT1G58340':'ABS4','AT5G52310':'COR78',
        'AT4G11290':'PRX39','AT5G09220':'AAP2','AT1G73480':'MAGL4'
    }
}

EXPANSION = [
 {'study':'GLDS-208','role':'direct root validation','assay':'microarray + RNA-seq','material':'root apex','priority':1,
  'rationale':'Closest independent tissue match to the GLDS-7 root discovery; tests tissue transfer across platforms.'},
 {'study':'GLDS-120','role':'context and genotype validation','assay':'RNA-seq','material':'multiple ecotypes, genotypes and light treatments','priority':2,
  'rationale':'Separates conserved response from ecotype, genotype and illumination effects.'},
 {'study':'GLDS-251','role':'gravity dose-response','assay':'RNA-seq','material':'seedlings under fractional gravity and blue light','priority':3,
  'rationale':'Supports a gravity-response curve rather than a binary flight-versus-ground interpretation.'},
 {'study':'GLDS-218','role':'regulatory mechanism','assay':'RNA-seq','material':'seedling development','priority':4,
  'rationale':'Adds alternative splicing as a regulatory layer beyond gene abundance.'},
 {'study':'GLDS-205','role':'mechanistic perturbation','assay':'microarray','material':'undifferentiated cells; HSFA2','priority':5,
  'rationale':'Tests heat-shock/redox regulation through an explicit transcription-factor perturbation.'},
 {'study':'GLDS-147','role':'gravity-sensing mechanism','assay':'microarray','material':'undifferentiated cells; ARG1','priority':6,
  'rationale':'Tests genetic involvement of a gravity-response component.'},
 {'study':'GLDS-213','role':'cell-autonomous transfer','assay':'microarray','material':'cell cultures flown on Shenzhou 8','priority':7,
  'rationale':'Tests whether part of the signature is retained without differentiated tissues.'},
 {'study':'GLDS-46','role':'radiation control','assay':'microarray','material':'seedlings exposed to gamma and HZE radiation','priority':8,
  'rationale':'Separates radiation-responsive genes from gravity-associated genes.'},
 {'study':'GLDS-136','role':'hypobaria control','assay':'microarray','material':'Arabidopsis under components of low pressure','priority':9,
  'rationale':'Separates low-pressure/hypoxia-like transcription from orbital exposure.'}
]

def main():
    bioc = pd.read_csv(SITE/'bioconductor-limma-results.csv')
    root = bioc[bioc.tissue.eq('root')].set_index('locus')
    external = external_sets()
    rows=[]
    genes=[]
    for panel, mapping in PANELS.items():
        for locus,symbol in mapping.items():
            g={'panel':panel,'locus':locus,'symbol':symbol,
               'glds7_root_log2fc':float(root.loc[locus,'log2fc']),
               'glds7_root_fdr':float(root.loc[locus,'fdr'])}
            genes.append(g)
        # Reference direction is the independently reproduced GLDS-7 root estimate.
        for study, frame in external.items():
            available=[g for g in mapping if g in frame.index]
            ref=np.array([np.sign(root.loc[g,'log2fc']) for g in available])
            obs=np.array([np.sign(frame.loc[g,'log2fc']) for g in available])
            same=int(np.sum(ref==obs)); n=len(available)
            p=stats.binomtest(same,n,.5,alternative='greater').pvalue if n else 1.0
            rows.append({'panel':panel,'study':study,'genes_available':n,'same_direction':same,
                         'direction_agreement':same/n if n else None,'exact_binomial_p':p,
                         'mean_external_log2fc':float(frame.loc[available,'log2fc'].mean()) if n else None,
                         'genes':','.join(f'{mapping[g]}:{frame.loc[g,"log2fc"]:.3f}' for g in available)})
    out=pd.DataFrame(rows)
    out.to_csv(SITE/'expanded-marker-validation.csv',index=False)
    pd.DataFrame(genes).to_csv(SITE/'expanded-marker-panels.csv',index=False)
    pd.DataFrame(EXPANSION).to_csv(SITE/'expanded-study-roadmap.csv',index=False)
    payload={'claim':'The GLDS-7 root panel is a tissue-specific discovery panel; the six-gene NASA panel is a literature-defined cross-study comparator, not a replacement.',
             'panels':{k:[{'locus':g,'symbol':s} for g,s in v.items()] for k,v in PANELS.items()},
             'existing_external_tests':rows,'expansion_studies':EXPANSION,
             'interpretation_boundary':'New-study roles are grounded in published study designs. Numeric validation reported here uses only the five external common-pipeline datasets already available in Barker et al. Supplementary Data 4.'}
    (SITE/'expanded-evidence.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8')
    print(out.to_string(index=False))

if __name__=='__main__': main()
