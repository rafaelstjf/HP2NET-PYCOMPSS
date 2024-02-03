from pycompss.api.task import task
from pycompss.api.binary import binary
from pycompss.api.parameter import DIRECTORY_INOUT

@binary(binary="raxmlHPC",
        working_dir=".",
        args= "-s {{alignment}} -m {{evo_model}} --HKY85 -f a -x {{seed}} -p {{seed}} -N {{bs_value}} -n {{out_suffix}} -w {{working_dir}}") # decorador para binários externos. O nome do programa deve ser passado como parâmetro
@task() # mesmo com o decorador tem-se que colocar o decorador de task
def raxml(alignment, evo_model, bs_value, out_suffix, seed, working_dir):
    pass

@task(basedir=DIRECTORY_INOUT)
def setup_phylip_data(basedir, test):
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
    import os, glob, tarfile, tarfile, shutil
    from Bio import AlignIO
    from pathlib import Path
    input_dir = os.path.join(basedir, 'input')
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
        pass # TODO raise an exception
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
                AlignIO.convert(f, "nexus", os.path.join(input_phylip_dir, f'{out_name}.fasta'), "fasta", molecule_type = "DNA")
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
       pass # TODO raise an exception
    return

