import subprocess
import asyncio
from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport

def query_gnomAD_local(chrom, pos, ref, alt, vcf_file):
    # Build bcftools query command
    query_cmd = f"""
    bcftools query -f '%CHROM %POS %REF %ALT %INFO/controls_popmax %INFO/controls_AC_popmax %INFO/controls_AN_popmax %INFO/controls_AF_popmax\n' {vcf_file}|
    awk '$1=="{chrom}" && $2=={pos} && $3=="{ref}" && $4=="{alt}"'
    """

    try:
        # Execute commands and capture output
        result = subprocess.run(query_cmd, shell=True, check=True, capture_output=True, text=True)
        output = result.stdout.strip()
        
        # Check if the output is empty
        if output:
            # Extract the value of Controls_popmax
            _, _, _, _, controls_popmax, controls_AC_popmax, controls_AN_popmax, controls_AF_popmax = output.split()
            return controls_popmax, controls_AC_popmax, controls_AN_popmax, controls_AF_popmax
        else:
            return None
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools query: {e}")
        return None

async def query_gnomAD_api(gene, chr, pos, ref, alt, cDNA):
    
    Max_MAF = -1
    Allele_Count = -1
    Allele_num = -1
    Max_AC = -1
    Max_AN = -1
    Continent = ""
    transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
    client = Client(transport=transport, fetch_schema_from_transport=True)
    
    # Define GraphQL queries
    query = gql(
        """
        query VariantsAFInGene($gene_symbol: String!) {
          gene(gene_symbol: $gene_symbol, reference_genome: GRCh37) {
            variants(dataset: gnomad_r2_1_controls) {
              chrom
              pos
              ref
              alt
              hgvsc
              genome {
                ac
                an
                af
                populations {
                  id
                  ac
                  an
                }
              }
            }
          }
        }
        """
    )
    variables = {"gene_symbol": gene}
    result = await client.execute_async(query, variable_values=variables)

    for variant in result["gene"]["variants"]:
        if variant['chrom'] == str(chr) and str(variant['pos']) == pos and variant['ref']== ref and variant['alt'] == alt and variant['hgvsc'] ==  cDNA:
            if variant['genome'] is not None:
              Allele_Count = variant['genome']['ac']
              Allele_num = variant['genome']['an']

              if variant['genome']['populations'] is not None:
                for var in variant['genome']['populations']:
                  if var['id'] == "afr":
                    # Africa
                    AFR_count = var['ac']
                    AFR_num = var['an']
                    AFR_AF = AFR_count/AFR_num if AFR_num != 0 else None
                  if var['id'] == "eas":
                    # East Asia
                    EAS_count = var['ac']
                    EAS_num = var['an']
                    EAS_AF = EAS_count/EAS_num if EAS_num != 0 else None
                  if var['id'] == "nfe":
                    # Non Finnish Europe
                    NFE_count = var['ac']
                    NFE_num = var['an']
                    NFE_AF = NFE_count/NFE_num if NFE_num != 0 else None
                  if var['id'] == "amr":
                    # Latino
                    AMR_count = var['ac']
                    AMR_num = var['an']
                    AMR_AF = AMR_count/AMR_num if AMR_num !=0 else None
                  if var['id'] == "sas":
                    # South Asia
                    SAS_count = var['ac']
                    SAS_num = var['an']
                    SAS_AF = var['ac']/var['an'] if SAS_num != 0 else None

                    Max_MAF = max(
                    -1 if AFR_AF is None else AFR_AF,
                    -1 if EAS_AF is None else EAS_AF,
                    -1 if NFE_AF is None else NFE_AF,
                    -1 if AMR_AF is None else AMR_AF,
                    -1 if SAS_AF is None else SAS_AF
                    )

                    # Identify the population with the highest AF
                    if Max_MAF == AFR_AF:
                      Max_AC = AFR_count
                      Max_AN = AFR_num
                      Continent = "African"
                    elif Max_MAF == EAS_AF:
                      Max_AC = EAS_count
                      Max_AN = EAS_num
                      Continent = "East Asian"
                    elif Max_MAF == NFE_AF:
                      Max_AC = NFE_count
                      Max_AN = NFE_num
                      Continent = "European (non-Finnish)"        
                    elif Max_MAF == AMR_AF:
                      Max_AC = AMR_count
                      Max_AN = AMR_num
                      Continent = "Latino"  
                    elif Max_MAF == SAS_AF:
                      Max_AC = SAS_count
                      Max_AN = SAS_num
                      Continent = "South Asian"  
                    elif Max_MAF == -1:
                      Max_AC = -1
                      Max_AC = -1   

    return Max_MAF, Allele_Count, Allele_num, Max_AC, Max_AN, Continent
