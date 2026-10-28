import os

RUNTIME = os.environ.get("HP2NET_RUNTIME", "parsl").lower()

if RUNTIME == "pycompss":
    from pycompss.api.api import compss_barrier, compss_wait_on

    def wait_for_all(futures):
        return compss_wait_on(futures)
else:
    import parsl
    from infra_manager import workflow_config, wait_for_all

import apps
import glob
import bioconfig
import os
import logging
import argparse
import math
from utils import CircularList

reuse = False
cache = dict()
# LOGGING SECTION
logger = logging.getLogger()
logging.basicConfig(level=logging.CRITICAL)


def build_raxml_trees(bio_config, basedir, preparation):
    max_workers = bio_config.workflow_core*bio_config.workflow_node
    ret_tree = list()
    pool = CircularList(math.floor(
        max_workers/int(bio_config.raxml_threads)))
    dir_ = os.path.join(os.path.join(basedir['dir'], "input"), "phylip")
    for input_file in glob.glob(os.path.join(dir_, '*.phy')):
        ret = apps.raxml(basedir=basedir,
                         config=bio_config,
                         inputs=preparation,
                         input_file=input_file, next_pipe=pool.next())
        pool.current(ret)
        ret_tree.append(ret)
    return apps.setup_tree_output(
        basedir=basedir, config=bio_config, inputs=ret_tree)


def raxml_snaq(bio_config, basedir, preparation, tree_output=None):
    max_workers = bio_config.workflow_core*bio_config.workflow_node
    result = list()
    ret_sad = (tree_output if tree_output is not None
               else build_raxml_trees(bio_config, basedir, preparation))
    logging.info("Using the Maximum Pseudo Likelihood Method")
    ret_ast = apps.astral(basedir, config=bio_config, inputs=[ret_sad])
    pool_phylo = CircularList(math.floor(
        max_workers/int(bio_config.snaq_threads)))
    for h in bio_config.snaq_hmax:
        ret_snq = apps.snaq(basedir, config=bio_config, hmax=h, inputs=[
                            ret_ast], next_pipe=pool_phylo.next())
        pool_phylo.current(ret_snq)
        result.append(ret_snq)
    return result


def raxml_phylonet(bio_config, basedir, preparation, tree_output=None):
    result = list()
    max_workers = bio_config.workflow_core*bio_config.workflow_node
    ret_sad = (tree_output if tree_output is not None
               else build_raxml_trees(bio_config, basedir, preparation))
    ret_rooted = apps.root_tree(basedir, config=bio_config, inputs=[ret_sad])
    logging.info("Using the Maximum Parsimony Method")
    out_dir = os.path.join(basedir['dir'], bio_config.phylonet_dir)
    pool_phylo = CircularList(math.floor(
        max_workers/int(bio_config.phylonet_threads)))
    for h in bio_config.phylonet_hmax:
        ret_spd = apps.setup_phylonet_data(
            basedir, config=bio_config, hmax=h, inputs=[ret_rooted])
        filename = os.path.join(
            out_dir, (basedir['tree_method'] + '_' + h + '_' + bio_config.phylonet_input))
        ret_phylonet = apps.phylonet(basedir, config=bio_config, input_file=filename, inputs=[
                                     ret_spd], next_pipe=pool_phylo.next())
        pool_phylo.current(ret_phylonet)
        result.append(ret_phylonet)
    return result


def build_iqtree_trees(bio_config, basedir, preparation):
    max_workers = bio_config.workflow_core*bio_config.workflow_node
    ret_tree = list()
    pool = CircularList(math.floor(
        max_workers/int(bio_config.iqtree_threads)))
    dir_ = os.path.join(os.path.join(basedir['dir'], "input"), "phylip")
    for input_file in glob.glob(os.path.join(dir_, '*.phy')):
        ret = apps.iqtree(basedir=basedir, config=bio_config,
                          inputs=preparation,
                          input_file=input_file, next_pipe=pool.next())
        pool.current(ret)
        ret_tree.append(ret)
    return apps.setup_tree_output(
        basedir=basedir, config=bio_config, inputs=ret_tree)


def iqtree_snaq(bio_config, basedir, preparation, tree_output=None):
    max_workers = bio_config.workflow_core*bio_config.workflow_node
    result = list()
    ret_sad = (tree_output if tree_output is not None
               else build_iqtree_trees(bio_config, basedir, preparation))
    logging.info("Using the Maximum Pseudo Likelihood Method")
    ret_ast = apps.astral(basedir, bio_config, inputs=[ret_sad])
    pool_phylo = CircularList(math.floor(
        max_workers/int(bio_config.snaq_threads)))
    for h in bio_config.snaq_hmax:
        ret_snq = apps.snaq(basedir, bio_config, h, inputs=[
                            ret_ast], next_pipe=pool_phylo.next())
        pool_phylo.current(ret_snq)
        result.append(ret_snq)
    return result


def iqtree_phylonet(bio_config, basedir, preparation, tree_output=None):
    result = list()
    max_workers = bio_config.workflow_core*bio_config.workflow_node
    ret_sad = (tree_output if tree_output is not None
               else build_iqtree_trees(bio_config, basedir, preparation))
    ret_rooted = apps.root_tree(basedir, bio_config, inputs=[ret_sad])
    logging.info("Using the Maximum Parsimony Method")
    out_dir = os.path.join(basedir['dir'], bio_config.phylonet_dir)
    pool_phylo = CircularList(math.floor(
        max_workers/int(bio_config.phylonet_threads)))
    for h in bio_config.phylonet_hmax:
        ret_spd = apps.setup_phylonet_data(
            basedir, bio_config, h, inputs=[ret_rooted])
        filename = os.path.join(
            out_dir, (basedir['tree_method'] + '_' + h + '_' + bio_config.phylonet_input))
        ret_phylonet = apps.phylonet(basedir, bio_config, filename, inputs=[
                                     ret_spd], next_pipe=pool_phylo.next())
        pool_phylo.current(ret_phylonet)
        result.append(ret_phylonet)
    return result


def mrbayes_snaq(bio_config, basedir, prepare_to_run):
    max_workers = bio_config.workflow_core*bio_config.workflow_node
    result = list()
    ret_tree = list()
    datalist = list()
    # append the input files
    dir_ = os.path.join(os.path.join(basedir['dir'], "input"), "nexus")
    datalist = glob.glob(os.path.join(dir_, '*.nex'))
    ret_mbsum = list()
    for input_file in datalist:
        ret_mb = apps.mrbayes(basedir, bio_config,
                              input_file=input_file, inputs=prepare_to_run)
        ret_mbsum.append(apps.mbsum(basedir, bio_config,
                         input_file=input_file, inputs=[ret_mb]))
    ret_pre_bucky = apps.setup_bucky_data(
        basedir, bio_config, inputs=ret_mbsum)
    wait_for_all([ret_pre_bucky])
    bucky_folder = os.path.join(basedir['dir'], "bucky")
    prune_trees = glob.glob(os.path.join(bucky_folder, "*.txt"))
    ret_bucky = list()
    for prune_tree in prune_trees:
        ret_bucky.append(apps.bucky(basedir, bio_config,
                         prune_file=prune_tree, inputs=[ret_pre_bucky]))
    ret_post_bucky = apps.setup_bucky_output(
        basedir, bio_config, inputs=ret_bucky)
    ret_pre_qmc = apps.setup_qmc_data(
        basedir, bio_config, inputs=[ret_post_bucky])
    ret_qmc = apps.quartet_maxcut(basedir, bio_config, inputs=[ret_pre_qmc])
    ret_tree.append(apps.setup_qmc_output(
        basedir, bio_config, inputs=[ret_qmc]))
    logging.info("Using the Maximum Pseudo Likelihood Method")
    pool_phylo = CircularList(math.floor(
        max_workers/int(bio_config.snaq_threads)))
    for h in bio_config.snaq_hmax:
        ret_snq = apps.snaq(basedir, bio_config, h,
                            inputs=ret_tree, next_pipe=pool_phylo.next())
        pool_phylo.current(ret_snq)
        result.append(ret_snq)
    return result

def _required_folders(config, basedir):
    """Return the output folders needed by one configured pipeline."""
    network_method = basedir['network_method']
    tree_method = basedir['tree_method']
    if network_method == 'MPL':
        if tree_method == 'RAXML':
            return [config.raxml_dir, config.astral_dir, config.snaq_dir]
        if tree_method == 'IQTREE':
            return [config.iqtree_dir, config.astral_dir, config.snaq_dir]
        if tree_method == 'MRBAYES':
            return [config.mrbayes_dir, config.bucky_dir, config.mbsum_dir,
                    config.quartet_maxcut_dir, config.snaq_dir]
    elif network_method == 'MP':
        if tree_method == 'RAXML':
            return [config.raxml_dir, config.phylonet_dir]
        if tree_method == 'IQTREE':
            return [config.iqtree_dir, config.phylonet_dir]
    return []


def prepare_to_run(config):
    """Prepare each dataset and its union of output folders exactly once."""
    preparations = list()
    datasets = dict()
    for basedir in config.workload:
        dataset = datasets.setdefault(
            basedir['dir'], {'basedir': basedir, 'folders': set()})
        dataset['folders'].update(_required_folders(config, basedir))

    for dataset in datasets.values():
        basedir = dataset['basedir']
        preparations.append(apps.setup_phylip_data(basedir, config))
        preparations.append(apps.create_folders(
            basedir, config, folders=sorted(dataset['folders'])))
    return preparations


def group_workloads(workload):
    """Group branches that can share the same tree-estimation futures."""
    groups = dict()
    for basedir in workload:
        key = (os.path.normpath(basedir['dir']), basedir['tree_method'])
        group = groups.setdefault(key, dict())
        group.setdefault(basedir['network_method'], basedir)
    return groups.values()


def main(**kwargs):
    logging.info('Starting the Workflow Orchestration')
    if kwargs["config_file"] is not None:
        config_file = kwargs["config_file"]
    else:
        config_file = 'default.ini'
    if kwargs["workload_file"] is not None:
        cf = bioconfig.ConfigFactory(
            config_file, custom_workload=kwargs["workload_file"])
    else:
        cf = bioconfig.ConfigFactory(config_file)
    bio_config = cf.build_config()
    logging.info(f"{hash(bio_config)}")
    if RUNTIME == "parsl":
        dkf_config = workflow_config(bio_config, **kwargs)
        parsl.load(dkf_config)
    results = list()
    prep = prepare_to_run(bio_config)

    if RUNTIME == "pycompss":
        compss_wait_on(prep)
    for branches in group_workloads(bio_config.workload):
        basedir = next(iter(branches.values()))
        tree_method = basedir['tree_method']

        if tree_method == 'RAXML':
            tree_output = build_raxml_trees(bio_config, basedir, prep)
            if 'MP' in branches:
                results.extend(raxml_phylonet(
                    bio_config, branches['MP'], prep, tree_output))
            if 'MPL' in branches:
                results.extend(raxml_snaq(
                    bio_config, branches['MPL'], prep, tree_output))
        elif tree_method == 'IQTREE':
            tree_output = build_iqtree_trees(bio_config, basedir, prep)
            if 'MP' in branches:
                results.extend(iqtree_phylonet(
                    bio_config, branches['MP'], prep, tree_output))
            if 'MPL' in branches:
                results.extend(iqtree_snaq(
                    bio_config, branches['MPL'], prep, tree_output))
        elif tree_method == 'MRBAYES':
            if 'MPL' in branches:
                results.extend(mrbayes_snaq(
                    bio_config, branches['MPL'], prep))
            if 'MP' in branches:
                logging.error('MRBAYES|MP is not a supported workflow')
        else:
            logging.error(f'Invalid tree method: {tree_method}')
    if bio_config.plot_networks == True:
        plot = apps.plot_networks(bio_config, inputs=results)
        wait_for_all([plot])
    else:
        wait_for_all(results)
    if RUNTIME == "pycompss":
        compss_barrier()
    else:
        parsl.dfk().cleanup()
    return


def cli():
    parser = argparse.ArgumentParser(
        description='Process the configuration file.')
    parser.add_argument('-s', '--settings', help='Settings file',
                        required=False, type=str, default=None)
    parser.add_argument('-w', '--workload', help='Workload file',
                        required=False, type=str, default=None)
    parser.add_argument(
        "-r", "--runinfo", help="Folder to store the Parsl logs", required=False, type=str, default=None)
    parser.add_argument(
        '-m', "--maxworkers", help="Max workers", required=False, type=int, default=None)

    args = parser.parse_args()

    main(config_file=args.settings, workload_file=args.workload,
         max_workers=args.maxworkers, runinfo=args.runinfo)


if __name__ == "__main__":
    cli()
