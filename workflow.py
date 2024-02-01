from pycompss.api.api import compss_wait_on
import apps

if __name__ == "__main__":
    r = apps.raxml(alignment="msa_artigo_denv.phy", evo_model="GTRGAMMA --HKY85",
                   bs_value="100", out_suffix="saida", seed="123")
    compss_wait_on(r)
