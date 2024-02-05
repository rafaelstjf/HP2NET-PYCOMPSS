from pycompss.api.api import compss_wait_on
import apps
import os
import bioconfig


def main():
    config_file = 'default.ini'
    cf = bioconfig.ConfigFactory(config_file)
    config = cf.build_config()
    basedir = (config.workload)[0]
    folder_list = [config.raxml_dir, config.astral_dir, config.snaq_dir]
    r = apps.create_folders(basedir, config, folder_list)
    compss_wait_on(r)
    genes = dict()
    genes = apps.convert_sequences(basedir, config)
    genes = compss_wait_on(genes)
    r_raxml = list()
    for g in genes["phylip"]:
        suffix = os.path.splitext(os.path.basename(g))[0]
        r_raxml.append(apps.raxml(alignment=g, evo_model=config.raxml_model,
                       bs_value=config.bootstrap, out_suffix=suffix, seed='123', working_dir=os.path.join(basedir["dir"], config.raxml_dir)))
    """
    r_sad = apps.setup_tree_output(basedir, config, r_raxml)
    r_pastral = apps.setup_astral(basedir, config, r_sad)
    r_astral = apps.astral(r_pastral["tree_output"], r_pastral["bs_file"], config.bootstrap, r_pastral["astral_output"])
    raxml_dir = os.path.join(basedir['dir'], config.raxml_dir)
    besttree_file = os.path.join(raxml_dir, config.raxml_output)
    spec_tree = os.path.join(basedir['dir'], os.path.join(config.astral_dir, config.raxml_dir))
    spec_tree = os.path.join(spec_tree, config.astral_output)
    output_folder = os.path.join(basedir['dir'], config.snaq_dir)
    r_snaq = apps.snaq(basedir["tree_method"], besttree_file, spec_tree, output_folder, "10", config.snaq_hmax, config.snaq_runs, r_astral)
    compss_wait_on(r_snaq)
    """


if __name__ == "__main__":
    main()
