from pycompss.api.api import compss_wait_on
import apps
import glob, os
def raxml_snaq(dataset):
    #get the list of genes
    input_folder = os.path.join(os.path.join(dataset,"input"), "phylip")
    raxml_folder = os.path.join(dataset, "raxml")
    genes = glob.glob(os.path.join(input_folder, "*.phy"))
    r_raxml = list()
    for gene in genes:
        suffix = os.path.splitext(os.path.basename(gene))[0]
        r = apps.raxml(alignment=gene, evo_model="GTRGAMMA",
                   bs_value="100", out_suffix=suffix, seed="123", working_dir=raxml_folder)
        r_raxml.append(r)
    apps.setup_phylip_data()
if __name__ == "__main__":
    raxml_snaq(os.path.join(os.getcwd(), "Denv_1_outgroup"))
    
    
