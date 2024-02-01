from pycompss.api.task import task
from pycompss.api.binary import binary
from pycompss.api.parameter import *

@binary(binary="raxmlHPC",
        working_dir=".",
        args= "-s {{alignment}} -m {{evo_model}} -f a -x {{seed}} -p {{seed}} -N {{bs_value}} -n {{out_suffix}}") # decorador para binários externos. O nome do programa deve ser passado como parâmetro
@task() # mesmo com o decorador tem-se que colocar o decorador de task
def raxml(alignment, evo_model, bs_value, out_suffix, seed):
    pass
