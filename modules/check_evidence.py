import numpy as np
from modules.autopvs1 import AutoPVS1, AutoPVS1CNV

# Splicing effect
def AffectSplicing(SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score):
    if SpliceAI_pred and ((dbscSNV_ada_score != '.' and float(dbscSNV_ada_score) >= 0.6) or (dbscSNV_rf_score != '.' and float(dbscSNV_rf_score) >= 0.6)):
        return True
    else:
        return False
def NotAffectSplicing(SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score):
    if not SpliceAI_pred and ((dbscSNV_ada_score != '.' and float(dbscSNV_ada_score) < 0.6) and (dbscSNV_rf_score != '.' and float(dbscSNV_rf_score) < 0.6)):
        return True
    else:
        return False

# Cancer hotspots: OS3、OM3、OP3
def isOS3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
    if isCancerHotspots:
        
        if CancerHotspots_sample >= 50 and CancerHotspots_count >= 10:
            return True
        else:
            return False
    elif not isCancerHotspots and isCOSMICHotspots:
        if COSMIC_sample >= 50 and COSMIC_count >= 10:
            return True
        else:
            return False
    elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
        if ExpertEvidenceCode == "OS3":
            return True
        else:
            return False
    else:
        return False


def isOM3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
    if isCancerHotspots:
        if CancerHotspots_sample < 50 and CancerHotspots_count >= 10:
            return True
        else:
            return False
    elif not isCancerHotspots and isCOSMICHotspots:
        if COSMIC_sample < 50 and COSMIC_count >= 10:
            return True
        else:
            return False
    elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
        if ExpertEvidenceCode == "OM3":
            return True
        else:
            return False
    else:
        return False

def isOP3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
    if isCancerHotspots:
        if CancerHotspots_count < 10:
            return True
        else:
            return False
    elif not isCancerHotspots and isCOSMICHotspots:
        if COSMIC_count < 10:
            return True
        else:
            return False
    elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
        if ExpertEvidenceCode == "OP3":
            return True
        else:
            return False
    else:
        return False

# Population data: SBVS1、SBS1、OP4 
def isSBVS1(inGnomAD, MAF, gene, Max_AC, Max_AN):
    if inGnomAD:
        if gene == "CDH1" and MAF > 0.002 and Max_AN > 2000 and Max_AC >= 5:
            return True
        elif gene == "PTEN" and MAF > 0.01 and Max_AN > 2000 and Max_AC >= 5:
            return True
        elif gene == "TP53" and MAF > 0.001 and Max_AC >= 5:
            return True
        elif gene not in ["CDH1", "PTEN", "TP53"] and MAF > 0.05:
            return True
        else:
            return False
    else:
        return False

def isSBS1(inGnomAD, MAF, gene, Max_AC, Max_AN):
    if inGnomAD:
        if gene == "CDH1" and (MAF <= 0.002 and MAF > 0.001) and Max_AN > 2000 and Max_AC >= 5:
            return True
        elif gene == "PTEN" and (MAF <= 0.01 and MAF > 0.001) and Max_AN > 2000 and Max_AC >= 5:
            return True
        elif gene == "TP53" and (MAF <= 0.001 and MAF > 0.0003) and Max_AC >= 5:
            return True
        elif gene not in ["CDH1", "PTEN", "TP53"] and (MAF <= 0.05 and MAF > 0.01):
            return True
        else:
            return False
    else:
        return False

def isOP4(inGnomAD, MAF, gene, Max_AC, Max_AN):
    
    if inGnomAD and gene == "CDH1":
        if Max_AC <= 2 and MAF <= 0.00001:
            return True
        elif Max_AC > 2 and MAF <= 0.00002:
            return True
        else:
            return False
    elif inGnomAD and gene == "PTEN":
        if Max_AC == 1 and MAF < 0.00001:
            return True
        elif Max_AC >= 2 and MAF <= 0.00002:
            return True
        else:
            return False
    elif gene == "TP53" and not inGnomAD:
        return True
    elif inGnomAD and gene not in ["CDH1", "PTEN", "TP53"]:
        if MAF < 0.0005:
            return True
        else:
            return False
    elif not inGnomAD and gene not in ["CDH1", "PTEN", "TP53"]:
        return True
    else:
        return False

# Computational Evidence: SBP1、OP1
def isSBP1(SIFT_pred,
          MutationAssessor_pred,
          FATHMM_pred,
          Polyphen2_HDIV_pred,
          Polyphen2_HVAR_pred,
          MutationTaster_pred,
          SpliceAI_pred,
          CADD_phred,
          REVEL_pred,
          dbscSNV_ada_score, dbscSNV_rf_score):
    # Evolutionary Conservation>=4, considered not to affect evolutionary conservation
    EvoConsPre_count = 0
    if SIFT_pred == 'T':
        EvoConsPre_count += 1
    if MutationAssessor_pred in ['L', 'N']:
        EvoConsPre_count += 1
    if FATHMM_pred == 'T':
        EvoConsPre_count += 1
    if Polyphen2_HDIV_pred == 'B':
        EvoConsPre_count += 1
    if Polyphen2_HVAR_pred == 'B':
        EvoConsPre_count += 1
    if MutationTaster_pred in ['N', 'P']:
        EvoConsPre_count += 1
    SPP1_count = 0
    if EvoConsPre_count >= 4:
        SPP1_count += 1
    if CADD_phred != '.' and float(CADD_phred) < 15:
        SPP1_count += 1
    if REVEL_pred != '.' and float(REVEL_pred) < 0.5:
        SPP1_count += 1
    if SPP1_count >= 2 or NotAffectSplicing(SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score):
    # At least two out of the three categories are considered to have no carcinogenic or splicing effects
        return True
    else:
        return False
    

def isOP1(SIFT_pred,
          MutationAssessor_pred,
          FATHMM_pred,
          Polyphen2_HDIV_pred,
          Polyphen2_HVAR_pred,
          MutationTaster_pred,
          SpliceAI_pred,
          CADD_phred,
          REVEL_pred,
          dbscSNV_ada_score, dbscSNV_rf_score):
    # Evolutionary Conservation>=4, considered to affect evolutionary conservation
    EvoConsPre_count = 0
    if SIFT_pred == 'D':
        EvoConsPre_count += 1
    if MutationAssessor_pred in ['H', 'M']:
        EvoConsPre_count += 1
    if FATHMM_pred == 'D':
        EvoConsPre_count += 1
    if Polyphen2_HDIV_pred in ['D', 'P']:
        EvoConsPre_count += 1
    if Polyphen2_HVAR_pred in ['D', 'P']:
        EvoConsPre_count += 1
    if MutationTaster_pred in ['A', 'D']:
        EvoConsPre_count += 1
    OP1_count = 0
    if EvoConsPre_count >= 4:
        OP1_count += 1
    if CADD_phred != '.' and float(CADD_phred) >= 15:
        OP1_count += 1
    if REVEL_pred != '.' and float(REVEL_pred) > 0.7:
        OP1_count += 1
    if OP1_count >= 2 or AffectSplicing(SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score): 
    # At least two out of the three categories are believed to have carcinogenic effects or splicing effects
        return True
    else:
        return False


    
# Functional evidence: OS2、SBS2
def isOS2(vcep_PS3, oncokb_O, oncokb_N, expertfunc_OS2, expertfunc_SBS2):
    # 1.Expert Interpretation
    if expertfunc_OS2:
        return True
    # 2.OncoKB
    elif not expertfunc_SBS2 and oncokb_O:
        return True
    # 3.ClinGen VCEP
    elif not expertfunc_SBS2 and not oncokb_N and vcep_PS3:
        return True
    else:
        return False

def isSBS2(vcep_BS3, oncokb_O, oncokb_N, expertfunc_OS2, expertfunc_SBS2):
    # 1.Expert Interpretation
    if expertfunc_SBS2:
        return True
    # 2.OncoKB
    elif not expertfunc_OS2 and oncokb_N:
        return True
    # 3.ClinGen VCEP
    elif not expertfunc_OS2 and not oncokb_O and vcep_BS3:
        return True
    else:
        return False
    

# Predictive data: SBP2、OVS1
def isSBP2(var_type, SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score):
    if var_type == "Synonymous" and (SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score):
        return True
    else:
        return False


def isOVS1(chr, pos, ref, alt, gDNA, isTSG):
    chr = str(chr)
    pos = str(pos)
    ref = str(ref)
    alt = str(alt)
    isOVS1 = False
    strength = ''
    consequence = ''
    if 'ins' not in gDNA and 'del' in gDNA:
        if '_' in gDNA:
            parts = pos.split('_')
            start = str(parts[0])
            end = str(parts[1])
            var = chr + '-' + start + '-' + end +'-DEL'
            demo = AutoPVS1CNV(var, 'hg19')
            consequence = demo.vep_consequence
            if demo.pvs1.strength_raw != "Unmet" and isTSG:
                isOVS1 = True
                strength = demo.pvs1.strength_raw
                return isOVS1, strength, consequence
            else:
                return isOVS1, strength, consequence
        else:# Single site deletion
            var = chr + '-' + pos + '-' + pos +'-DEL'
            demo = AutoPVS1CNV(var, 'hg19')
            consequence = demo.vep_consequence
            if demo.pvs1.strength_raw != "Unmet" and isTSG:
                isOVS1 = True
                strength = demo.pvs1.strength_raw
                return isOVS1, strength, consequence
            else:
                return isOVS1, strength, consequence
    elif 'dup' in gDNA:
        if '_' in gDNA:# Regional duplication
            parts = pos.split('_')
            start = str(parts[0])
            end = str(parts[1])
            var = chr + '-' + start + '-' + end +'-DUP'
            demo = AutoPVS1CNV(var, 'hg19')
            consequence = demo.vep_consequence
            if demo.pvs1.strength_raw != "Unmet" and isTSG:
                isOVS1 = True
                strength = demo.pvs1.strength_raw
                return isOVS1, strength, consequence
            else:
                return isOVS1, strength, consequence
        else:# Single site duplication
            var = chr + '-' + pos + '-' + pos +'-DUP'
            demo = AutoPVS1CNV(var, 'hg19')
            consequence = demo.vep_consequence
            if demo.pvs1.strength_raw != "Unmet" and isTSG:
                isOVS1 = True
                strength = demo.pvs1.strength_raw
                return isOVS1, strength, consequence
            else:
                return isOVS1, strength, consequence
    elif 'delins' not in gDNA and ref != '-' and alt != '-':
        var = chr + '-' + pos + '-' + ref + '-' + alt
        demo = AutoPVS1(var, 'hg19')
        consequence = demo.consequence
        if demo.islof and isTSG:
            isOVS1 = True
            strength = demo.pvs1.strength_raw
            return isOVS1, strength, consequence
        else:
            return isOVS1, strength, consequence
    else:
        return isOVS1, strength, consequence

# OM1、OM2
def isOM1(domain_name):
    return domain_name

def isOM2(isONG, isTSG, var_type):
    if (isONG or isTSG) and ("InFrameDeletion" in var_type or "InFrameInsertion" in var_type):
        return True
    elif isTSG and ("CdsStopDeletion" in var_type):
        return True
    else:
        return False