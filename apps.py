from pycompss.api.task import task
from pycompss.api.binary import binary
from pycompss.api.parameter import *
from bioconfig import BioConfig
import os, glob, tarfile, logging, tarfile, shutil
from pycompss.api.julia import julia
from Bio import AlignIO
from pathlib import Path
from appsexception import *
@binary(binary="raxmlHPC",
        working_dir=".",
        args= "-s {{alignment}} -m {{evo_model}} --HKY85 -f a -x {{seed}} -p {{seed}} -N {{bs_value}} -n {{out_suffix}} -w {{working_dir}}") # decorador para binários externos. O nome do programa deve ser passado como parâmetro
@task() # mesmo com o decorador tem-se que colocar o decorador de task
def raxml(alignment, evo_model, bs_value, out_suffix, seed, working_dir):
    pass

@task(basedir=IN, config=IN,returns=1)
def convert_sequences(basedir: dict, config: BioConfig):
    """Extract the sequence alignments tar file and convert the gene alignments from the nexus format to the phylip format.

    Parameters:
            basedir: it is going to search for a tar file with nexus files. The script will create:
            seqdir=input/nexus
            seqdir=input/phylip
    Returns:
        returns an parsl's AppFuture.

    Raises:
        PhylipMissingData --- if cannot find a tar file with nexus files.


    TODO: Provide provenance.
    NB:
        Stdout and Stderr are defaulted to parsl.AUTO_LOGNAME, so the log will be automatically 
        named according to task id and saved under task_logs in the run directory.
    """

    logging.info(f'Converting Nexus files to Phylip on {basedir["dir"]}')
    input_dir = os.path.join(basedir['dir'], 'input')
    sequence_dir = os.path.join(input_dir, 'sequence')
    input_format = 0
    # First the sequences are extracted
    Path(sequence_dir).mkdir(exist_ok=True)
    if os.path.exists(sequence_dir):
        tar_file = basedir['sequences']
        tar = tarfile.open(tar_file, "r:gz")
        tar.extractall(path=sequence_dir)
    # Now one file is opened to check its format
    sequences = glob.glob(os.path.join(sequence_dir, '*'))
    if len(sequences) == 0:
        raise AlignmentConversion(input_phylip_dir)
    with open(sequences[0], 'r') as s_file:
        line = s_file.readline()
        if "#NEXUS" in line:
            input_format = 0 # nexus
        elif ">" in line:
            input_format = 1 # fasta
        else:
            input_format = 2 # other .i.e. phylip

    input_nexus_dir = os.path.join(input_dir, 'nexus')
    input_phylip_dir = os.path.join(input_dir, 'phylip')
    input_fasta_dir = os.path.join(input_dir, 'fasta')
    # So, some work must be done. Build the Nexus directory
    if not os.path.isdir(input_nexus_dir):
        Path(input_nexus_dir).mkdir(exist_ok=True)
    if not os.path.isdir(input_phylip_dir):
        Path(input_phylip_dir).mkdir(exist_ok=True)
    if not os.path.isdir(input_fasta_dir):
        Path(input_fasta_dir).mkdir(exist_ok=True)
    # Now, use the function to convert nexus to phylip.
    files = glob.glob(os.path.join(sequence_dir,'*'))
    try:
        for f in files:
            out_name = os.path.basename(f).split('.')[0]
            if input_format == 0:
                AlignIO.convert(f, "nexus", os.path.join(input_phylip_dir, f'{out_name}.phy'), "phylip-sequential", molecule_type = "DNA")
                AlignIO.convert(f, "nexus", os.path.join(input_fasta_dir, f'{out_name}.fasta'), "fasta", molecule_type = "DNA")
                shutil.copyfile(f, os.path.join(input_nexus_dir, os.path.basename(f)))
            if input_format == 1:
                AlignIO.convert(f, "fasta", os.path.join(input_phylip_dir, f'{out_name}.phy'), "phylip-sequential", molecule_type = "DNA")
                AlignIO.convert(f, "fasta", os.path.join(input_nexus_dir, f'{out_name}.nex'), "nexus", molecule_type = "DNA")
                shutil.copyfile(f, os.path.join(input_fasta_dir, os.path.basename(f)))
            if input_format == 2:
                AlignIO.convert(f, "phylip-sequential", os.path.join(input_nexus_dir, f'{out_name}.nex'), "nexus", molecule_type = "DNA")
                AlignIO.convert(f, "phylip-sequential", os.path.join(input_fasta_dir, f'{out_name}.fasta'), "fasta", molecule_type = "DNA")
                shutil.copyfile(f, os.path.join(input_phylip_dir, os.path.basename(f)))
    except Exception as e:
        raise AlignmentConversion(input_phylip_dir)
    seq_dict = dict()
    seq_dict["nexus"] = list()
    seq_dict["fasta"] = list()
    seq_dict["phylip"] = list()
    for f in files:
            out_name = os.path.basename(f).split('.')[0]
            seq_dict["nexus"].append(os.path.join(input_nexus_dir, f'{out_name}.nex'))
            seq_dict["fasta"].append(os.path.join(input_fasta_dir, f'{out_name}.fasta'))
            seq_dict["phylip"].append(os.path.join(input_phylip_dir, f'{out_name}.phy'))
    return seq_dict

@task(basedir=IN, config=IN, inputs=COLLECTION_INOUT)
def setup_tree_output(basedir: dict,
                      config: BioConfig,
                      inputs: dict
                      ):
    """Create the phylogenetic tree software (raxml, iqtree,...) best tree file and organize the temporary files to subsequent softwares 

    Parameters:
        basedir: current working directory
    Returns:
        returns an parsl's AppFuture

    TODO: 
        Provide provenance.

    NB:
        Stdout and Stderr are defaulted to parsl.AUTO_LOGNAME, so the log will be automatically 
        named according to task id and saved under task_logs in the run directory.
    """
    work_dir = basedir['dir']
    tree_method = basedir['tree_method']
    logging.info(f'Setting up the tree output on {work_dir}')
    if(tree_method == "RAXML"):
        raxml_dir = os.path.join(work_dir, config.raxml_dir)
        bootstrap_dir = os.path.join(raxml_dir, "bootstrap")
        besttree_file = os.path.join(raxml_dir, config.raxml_output)
        try:
            Path(bootstrap_dir).mkdir(exist_ok=True)
        except Exception:
            raise FolderCreationError(bootstrap_dir)
        old_files = glob.glob(f'{bootstrap_dir}/*')
        try:
            for f in old_files:
                os.remove(f)
        except Exception:
            raise FolderDeletionError(bootstrap_dir)
        try:
            files = glob.glob(os.path.join(raxml_dir,'RAxML_bootstrap.*'))
            for f in files:
                os.rename(f, os.path.join(bootstrap_dir, os.path.basename(f)))
            # compress and remove the bootstrap files
            with tarfile.open(os.path.join(raxml_dir, "contrees.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(raxml_dir,'RAxML_bipartitions.*'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            raxml_input = open(besttree_file, 'w')
            files = glob.glob(os.path.join(raxml_dir, 'RAxML_bestTree.*'))
            trees = ""
            for f in files:
                gen_tree = open(f, 'r')
                trees += gen_tree.readline()
                gen_tree.close()
            raxml_input.write(trees)
            raxml_input.close()
            with tarfile.open(os.path.join(raxml_dir, "besttrees.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(raxml_dir, 'RAxML_bestTree.*'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(raxml_dir, "bipartitionsBranchLabels.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(raxml_dir, 'RAxML_bipartitionsBranchLabels.*'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(raxml_dir, "info.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(raxml_dir, 'RAxML_info.*'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
        except IOError:
            raise FileCreationError(raxml_dir)
    elif(tree_method == "IQTREE"):
        phylip_dir = os.path.join(work_dir, os.path.join("input", "phylip"))
        iqtree_dir = os.path.join(work_dir, config.iqtree_dir)
        besttree_file = os.path.join(iqtree_dir, config.iqtree_output)
        try:
            files = glob.glob(os.path.join(phylip_dir, '*.iqtree'))
            files += glob.glob(os.path.join(phylip_dir, '*.treefile'))
            files += glob.glob(os.path.join(phylip_dir, '*.mldist'))
            files += glob.glob(os.path.join(phylip_dir, '*.nex'))
            files += glob.glob(os.path.join(phylip_dir, '*.contree'))
            files += glob.glob(os.path.join(phylip_dir, '*.log'))
            files += glob.glob(os.path.join(phylip_dir, '*.ckp.gz'))
            files += glob.glob(os.path.join(phylip_dir, '*.bionj'))
            #files += glob.glob(os.path.join(phylip_dir, '*.reduced'))
            #files += glob.glob(os.path.join(phylip_dir, '*.boottrees'))
            files += glob.glob(os.path.join(phylip_dir, '*.ufboot'))
            for f in files:
                new_f = os.path.join(iqtree_dir, os.path.basename(f))
                os.replace(f, new_f)
            iq_input = open(besttree_file, 'w+')
            files = glob.glob(os.path.join(iqtree_dir, '*.treefile'))
            trees = ""
            for f in files:
                gen_tree = open(f, 'r')
                trees += gen_tree.readline()
                gen_tree.close()
            iq_input.write(trees)
            iq_input.close()
        except IOError:
            raise FileCreationError(iqtree_dir)
        bootstrap_dir = os.path.join(iqtree_dir, "bootstrap")
        try:
            Path(bootstrap_dir).mkdir(exist_ok=True)
        except Exception:
            raise FolderCreationError(bootstrap_dir)
        old_files = glob.glob(f'{bootstrap_dir}/*')
        try:
            for f in old_files:
                os.remove(f)
        except Exception:
            raise FolderDeletionError(bootstrap_dir)
        try:
            with tarfile.open(os.path.join(iqtree_dir, "iqtree.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.iqtree'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "treefile.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.treefile'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "mldist.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.mldist'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "nex.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.nex'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "contree.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.contree'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "log.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.log'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "ckp.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.ckp.gz'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "bionj.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.bionj.gz'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            with tarfile.open(os.path.join(iqtree_dir, "reduced.tgz"), "w:gz") as tar:
                files = glob.glob(os.path.join(iqtree_dir,'*.reduced.gz'))
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
                for f in files:
                    os.remove(f)
            files = glob.glob(os.path.join(iqtree_dir,'*.ufboot'))
            for f in files:
                os.rename(f, os.path.join(bootstrap_dir, os.path.basename(f)))
        except:
            raise FileCreationError(iqtree_dir)
    return

@task(basedir=IN, config=IN, setup=IN)
def setup_astral(basedir, config, setup):

    work_dir = basedir['dir']
    tree_method = basedir['tree_method']
    mapping = basedir['mapping']
    logging.info(f'ASTRAL called with {work_dir}')
    astral_dir = os.path.join(work_dir,config.astral_dir)
    tree_output = ""
    astral_output = ""
    if(tree_method == "RAXML"):
        try:
            astral_raxml = os.path.join(astral_dir, config.raxml_dir)
            Path(astral_raxml).mkdir(exist_ok=True)
        except Exception:
            print("Failed to create the raxml bootstrap folder!")
        bs_file = os.path.join(astral_raxml,'BSlistfiles')
        raxm_dir = os.path.join(work_dir, config.raxml_dir)
        tree_output = os.path.join(raxm_dir,config.raxml_output)
        boot_strap = os.path.join(os.path.join(work_dir,config.raxml_dir),"bootstrap/*")
        with open(bs_file, 'w') as f:
            for i in glob.glob(boot_strap):
                f.write(f'{i}\n')
        astral_output = os.path.join(astral_raxml, config.astral_output)
    elif(tree_method == "IQTREE"):
        try:
            astral_iqtree = os.path.join(astral_dir, config.iqtree_dir)
            Path(astral_iqtree).mkdir(exist_ok=True)
        except Exception:
            print("Failed to create the raxml bootstrap folder!")
        bs_file = os.path.join(astral_iqtree,'BSlistfiles')
        iqtree_dir = os.path.join(work_dir, config.iqtree_dir)
        tree_output = os.path.join(iqtree_dir,config.iqtree_output)
        boot_strap = os.path.join(os.path.join(work_dir,config.iqtree_dir),"bootstrap/*")
        with open(bs_file, 'w') as f:
            for i in glob.glob(boot_strap):
                f.write(f'{i}\n')
        astral_output = os.path.join(astral_iqtree, config.astral_output)
    # Return to Parsl to be executed on the workflow
    params = f'-i {tree_output} -b {bs_file} -r {config.bootstrap} -o {astral_output}'
    if len(mapping) > 0:
        map_filename = os.path.join(astral_dir, 'mapping.dat')
        with open(map_filename, 'w') as map_:
            species = mapping.split(';')
            for specie in species:
                map_.write(specie.strip() + '\n')
            map_.close()
    r_dict = dict()
    r_dict["tree_output"] = tree_output
    r_dict["bs_file"] = bs_file
    r_dict["astral_output"] = astral_output
    return r_dict

@binary(binary="java",
         working_dir='.',
         args="-jar /prj/posgrad/rafaelst/Documentos/GitHub/ASTRAL/Astral/astral.5.7.8.jar -i {{tree_output}} -b {{bs_file}} -r {{bs_value}} -o {{astral_output}}")
@task()
def astral(tree_output, bs_file, bs_value, astral_output):
    pass

@binary(binary="julia",
        working_dir = '.',
        args="/prj/posgrad/rafaelst/Documentos/GitHub/HP2NET-PYCOMPSS/scripts/snaq.jl {{tree_method}} {{gen_tree}} {{spec_tree}} {{output_folder}} {{num_threads}} {{hmax}} {{runs}}")
@task(astral=IN)
def snaq(tree_method, gen_tree, spec_tree, output_folder, num_threads, hmax, runs, astral):
    pass

@task(basedir=IN, config=IN, folders=COLLECTION_IN)
def create_folders(basedir: dict,
                   config: BioConfig,
                   folders=[]):
    work_dir = basedir['dir']
    logging.info(f'Removing folders from old executions')
    for folder in folders:
        full_path = os.path.join(work_dir, folder)
        if(os.path.exists(full_path)):
            try:
                shutil.rmtree(full_path, ignore_errors=True)
            except Exception:
                raise FolderDeletionError(full_path)
    logging.info(f'Creating folders in {work_dir}')
    for folder in folders:
        full_path = os.path.join(work_dir, folder)
        try:
            Path(full_path).mkdir(exist_ok=True)
        except Exception:
            FolderCreationError(full_path)
    return
