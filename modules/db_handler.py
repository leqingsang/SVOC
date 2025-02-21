import pymysql, sqlite3, re
from warnings import filterwarnings

global db, cursor

db = sqlite3.connect("svocdb/SVOC_v20250108.db", uri=True)
cursor = db.cursor()



# Execute SQL command, capture errors and warnings
def pymysql_cursor(sql, params=None):
  # try... except... can only capture errors
  filterwarnings("error", category = pymysql.Warning)
  try:
    # cursor.execute(sql) will automatically return the number of affected rows
    if params is not None:
        effected_rows = cursor.execute(sql, params)
    else:
        effected_rows = cursor.execute(sql)
    db.commit()
  except pymysql.Error as e:
    db.rollback()
    raise e  
  except pymysql.Warning as w:
    db.rollback()
    raise w  
  else:
    if re.match('^SELECT', sql):
      # cursor.fetchone() return tuple; cursor.execute(sql) returns the number of rows affected 
      # one query result: ((1,)); multiple results: ((1,),(2,))
      result = cursor.fetchall()
      if len(result) == 1:
        single_result = result[0][0]
        return single_result
      elif len(result) > 1:
        multiple_results = []
        for ele in result:
          multiple_results.append(ele[0])
        return multiple_results
    elif re.match('^LOAD DATA LOCAL INFILE', sql):
      return cursor.rowcount

def getHotspots(gene, refAA, altAA, posAA, aa):

    isExist = False
    sample = 0
    count = 0
    sample = pymysql_cursor("SELECT Mutation_Count from SNVHotspots WHERE Gene = ? and Amino_Acid_Position = ? and Reference_Amino_Acid LIKE ? and Variant_Amino_Acid LIKE ?;", (gene, posAA, f"{refAA}%", f"{altAA}%"))
    raw_count = pymysql_cursor("SELECT Variant_Amino_Acid from SNVHotspots WHERE Gene = ? and Amino_Acid_Position = ? and Reference_Amino_Acid LIKE ? and Variant_Amino_Acid LIKE ?;", (gene, posAA, f"{refAA}%", f"{altAA}%"))
    isExist = bool(raw_count)
    if isExist:
      count = int(re.findall(r'\d+', raw_count)[0])
    else: # Search Indel
      
      sample = pymysql_cursor("SELECT Mutation_Count from INDELHotspots WHERE Gene = ? and Variant_Amino_Acid LIKE ?;", (gene, f"{aa}%"))
      raw_count = pymysql_cursor("SELECT Variant_Amino_Acid from INDELHotspots WHERE Gene = ? and Variant_Amino_Acid LIKE ?;", (gene, f"{aa}%"))
      isExist = bool(raw_count)
      if isExist:
        count = int(re.findall(r':(\d+)', raw_count)[0])

    return isExist, sample, count

def getCOSMICHotspots(gene, AAchange_single):
    refAA = re.search(r"p\.([A-Za-z]|\*)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(1)
    posAA = re.search(r"p\.([A-Za-z]|\*|)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(2)
    AAchange_position = 'p.'+refAA+posAA
    isExist = False
    sample = 0
    count = 0
    sample = pymysql_cursor("SELECT SUM(COSMIC_SAMPLE_MUTATED) from Mutation WHERE Gene = ? AND Mutation_AA LIKE ?;", (gene, f"{AAchange_position}_"))
    if '_' in AAchange_position: # When multiple amino acid changes occur, usually followed by a single or multiple characters, the wildcard '%' can be used
       sample = pymysql_cursor("SELECT SUM(COSMIC_SAMPLE_MUTATED) from Mutation WHERE Gene = ? AND Mutation_AA LIKE ?;", (gene, f"{AAchange_position}%"))
    count = pymysql_cursor("SELECT SUM(COSMIC_SAMPLE_MUTATED) from Mutation WHERE Gene = ? AND Mutation_AA = ?;", (gene, AAchange_single))
    isExist = bool(count)
    return isExist, sample, count

def getExpertCuratedHotspots(gene, gDNA):
   isExist = False
   ExpertEvidenceCode = ''
   ExpertEvidenceCode = pymysql_cursor("SELECT EvidenceCode from ExpertCuratedHotspots WHERE Gene = ? AND gDNA = ?;", (gene, gDNA))
   isExist = bool(ExpertEvidenceCode)
   ExpertCuratedHotspots_Basis = pymysql_cursor("SELECT Basis from ExpertCuratedHotspots WHERE Gene = ? AND gDNA = ?;", (gene, gDNA))

   return isExist, ExpertEvidenceCode, ExpertCuratedHotspots_Basis

def getVCEPMetCodes(Transcript, cDNA):
    hgvs = str(Transcript) + ":" + str(cDNA)
    MetCodes = pymysql_cursor("SELECT MetCodes from ClinGenVCEP WHERE ? LIKE HGVS;", (f"%{hgvs}%",))
    return MetCodes
def getOncoKBRes(aa,gene):
    oncokb_class = pymysql_cursor("SELECT Oncogenic_effect from OncoKBVariant where Alteration = ? and Gene = ?;", (aa,  gene))
    oncokb_Description = pymysql_cursor("SELECT Description from OncoKBVariant where Alteration = ? and Gene = ?;", (aa,  gene))
    return oncokb_class, oncokb_Description

def getExpertFuncRes(aa,gene,gDNA):
    if aa:
      expertfunc_class = pymysql_cursor("SELECT Oncogenic_effect from ExpertFuncRes where Alteration = ? and Gene = ?;", (aa,  gene))
      expertfunc_Description = pymysql_cursor("SELECT Description from ExpertFuncRes where Alteration = ? and Gene = ?;", (aa,  gene))
    else:
      expertfunc_class = pymysql_cursor("SELECT Oncogenic_effect from ExpertFuncRes where Alteration = ? and Gene = ?;", (gDNA,  gene))
      expertfunc_Description = pymysql_cursor("SELECT Description from ExpertFuncRes where Alteration = ? and Gene = ?;", (gDNA,  gene))
    return expertfunc_class, expertfunc_Description

def isOncoGene(gene):
    role_in_cancer = pymysql_cursor("SELECT RoleinCancer from Gene where Gene = ?;", (gene,))

    if role_in_cancer and "oncogene" in role_in_cancer:
       return True
    else:
       return False

def isTSG(gene):  
    role_in_cancer = pymysql_cursor("SELECT RoleinCancer from Gene where Gene = ?;", (gene,))
    if role_in_cancer and "TSG" in role_in_cancer:
       return True
    else:
       return False


def getOncokbDomain(gene,startAA,endAA):
   domain_name = ''
   domain_name = pymysql_cursor("SELECT Short_name from OncoKBDomain where Gene = ? and AAChange_start <= ? and AAChange_end >= ?;", (gene, startAA, endAA))
   return domain_name