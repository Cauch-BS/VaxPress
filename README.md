# VaxPress

![Build Status](https://github.com/Cauch-BS/VaxPress/actions/workflows/build.yaml/badge.svg)

VaxPress is a codon optimizer platform tailored for mRNA vaccine
development. It refines coding sequences starting from protein or
RNA sequences to boost both storage stability and *in vivo* protein
expression. Plus, additional properties can be easily programmed
into the optimization process with just a few lines of code via a
pluggable interface. For the detailed information about VaxPress,
including its options and algorithmic features, please refer to the
[VaxPress documentation page](https://vaxpress.readthedocs.io/).

VaxPress has been tested on Python 3.9, 3.10, 3.11 and 3.12. 

# Installation

### pip

You can install VaxPress via pip.

#### Installing

```bash
# Create a virtual environment for VaxPress
python -m venv /path/to/vaxpress-env

# Activate the virtual environment
source /path/to/vaxpress-env/bin/activate

# Install VaxPress alone
pip install vaxpress

# Alternatively, install VaxPress with LinearFold (only for non-commercial uses)
pip install 'vaxpress[nonfree]'
```

#### Running

```bash
# Activate the virtual environment
source /path/to/vaxpress-env/bin/activate

# Run VaxPress
vaxpress -h
```

#### iCodon Dependency

If you wish to activate the iCodon predicted stability
(`--iCodon-weight`) in the fitness function, ensure you have
working installations of *R,* *rpy2* (version >= 3.0) and
*iCodon.*  For detailed installation instructions, visit
[iCodon's GitHub page](https://github.com/santiago1234/iCodon/).

#### IDT Complexity Score

If you wish to use the gBlocks IDT complexity score, as accessible at
the [IDT site](https://sg.idtdna.com/site/order/gblockentry) use the argument
`--idt-complexity-weight` in the fitness function and export your personal
details as environmental variables. For example if you are using a UNIX based
operating system (such as Linux, MacOS) you can use the following commands:

```bash
export IDT_USER='your_IDT_username'
export IDT_PASSWD='your_IDT_password'
export IDT_CLIENT'your_IDT_client_id'
export IDT_API = 'your_IDT_client_secret'
```

where `'your_IDT_username'` should be replaced with your username. 

More information on how to use the IDT API is available at [SciTools API site](https://sg.idtdna.com/pages/tools/apidoc). Note the maximum number of allowed API requests per minute is 500, so it is generally not recommended that you use a population greater than 500 when using the IDT compleity score. 

For more information as to how to use the IDT complexity score, see the **Usage** section. 

### Conda

Alternatively, you may also install VaxPress via a conda package:

#### Installation

```bash
conda create -n vaxpress -y -c changlabsnu -c bioconda -c conda-forge vaxpress
```

#### Running

```bash
# Activate the environment
conda activate vaxpress

# Run VaxPress
vaxpress -h
```

### Singularity

To run VaxPress via Singularity, you will need to install the
[Singularity CE](https://sylabs.io/singularity/) first.
Download the container image from
[the GitHub project page](https://github.com/ChangLabSNU/VaxPress/releases)
and place it in a directory of your choice.

```bash
singularity run vaxpress.sif -h
```

When using the Singularity image, both the input and output must
be somewhere inside your home directory for VaxPress to run without
complicated directory binding configurations for Singularity.

# Usage

## Quick Start

Here's a basic command-line instruction to start using VaxPress.
Note that `-i` and `-o` options are mandatory:

```bash
vaxpress -i spike.fa -o output --iterations 1000 -p 32
```

### Input

VaxPress requires a FASTA format input file that contains the CDS
(CoDing Sequence) to be optimized. In case the FASTA file holds a
protein sequence, the additional `--protein` switch is required.

### Number of Iterations

The `--iterations` option is set to `10` by default. However,
for thorough optimization, it's recommended to use at least `500`
iterations. The optimal number of iterations may differ depending
on the length, composition of the input, and the selected optimization
settings. It's important to note that the optimization process may
stop before completing all the specified iterations if no progress
is observed over several consecutive cycles. Guidelines for setting
the appropriate number of iterations and other optimization parameters
can be found in the
[Tuning Optimization Parameters](https://vaxpress.readthedocs.io/en/latest/user_guides.html#tuning-parameters)
section.

You can set `--iterations` to `0` to generate VaxPress's sequence
evaluation report without any optimization.

### Multi-Core Support

You can use multiple CPU cores for optimization with the `-p` or
`--processes` option.

### More About Options

VaxPress offers the flexibility to adjust optimization strategies
in detail and integrate with LinearDesign. It also allows several
more convenient functions such as preset parameters, user-defined
custom scoring functions, and etc. For comprehensive explanation,
please refer to [the manual](https://vaxpress.readthedocs.io/en/latest/).

## Output

Once you've run VaxPress, the specified output directory will contain
the following five files:

- ``report.html``: A summary report detailing the result and
  optimization process.
- ``best-sequence.fasta``: The refined coding sequence.
- ``checkpoints.tsv``: The best sequences and the evaluation results
  at each iteration.
- ``log.txt``: Contains the logs that were displayed in the console.
- ``parameters.json``: Contains the parameters employed for the optimization.
  This file can be feeded to VaxPress with the `--preset` option to duplicate
  the set-up for other sequence.

## Using VaxPress Polish

To optimize VaxPress candidate sequences for production, the `vaxpress-polish` command can be used. There are two commands available for `vaxpress-polish`. The first is `find-complexity` and the second is `identify-alternatives`. Assuming you have activated your IDT account appropriately (for more information on this part, see the IDT complexity score section below insallation), you can use `find-complexity` to identify which parts of the sequence are increasing complexity. A complexity above 10 is considered too complex and will likely be costly to synthesize. Here is an example:

```bash
>>> vaxpress-polish find-complexity -i s901_temp.fasta

Complexity score:
 {
   'Overall Repeat': {   'Message': 'One or more repeated sequences greater than 8 bases comprise '
                                     '64% of the overall sequence. Solution: Redesign to reduce '
                                     'the repeats to be less than 40% of the sequence.',
                          'Score': 9.6}
}
```
Meanwhile the `identify-alternatives` identifies codons which can be changed with minimal change in VaxPress scores. The way it does this is by focusing on unpaired bases. The corresponding change in CAI due to the change in codon is written next to the option given. Here is an example:

```bash
>>> vaxpress-polish identify-alternatives -i s901_temp.fasta --cds-start 75 --cds-end 2243 --roi-start 1481 --roi-end 1521
ROI Sequence:  AAAGAGAGGCUCUGGGGAGGGCCGGGGAAGCCUGCUGACCU
ROI Structure: ........))))))))).))))))))))...)).)))...)
        1479 agG   -0.13
        1479 CgC   -0.61
        1479 CgG   -0.30
        1479 CgU   -1.55
           1482 aaA   -0.21
              1485 Cga   -1.05
              1485 CgC   -0.61
              1485 CgG   -0.30
              1485 CgU   -1.55
              1485 agG   -0.13
              1485 CgC   -0.61
              1485 CgG   -0.30
              1485 CgU   -1.55
                                      1509 UCA   -0.48
                                      1509 UCc   -0.19
                                      1509 UCG   -2.29
                                      1509 UCU   -0.22
                                      1509 UCA   -0.48
                                      1509 UCc   -0.19
                                      1509 UCG   -2.29
                                      1509 UCU   -0.22
                                      1509 agU   -0.49
                                      1509 UCA   -0.48
                                      1509 UCG   -2.29
                                      1509 UCU   -0.22
                                         1512 cuA   -2.28
                                         1512 cuC   -1.02
                                         1512 cuU   -1.36
                                         1512 UuA   -2.05
                                               1518 acA   -0.11
                                               1518 acG   -1.67
                                               1518 acU   -0.32
Results written to codon_polish_results.fasta
```

For the choice `1479 agG -0.13` this means that the codon starting from the 1479th base (`AGA`) can be exhanged fo `AGG` with a final CAI change of `-0.13`. 

The default output (also modifiable through the `-o` option) is `codon_polish_results.fasta` which contains the changed secondary stucture and the changed sequence given that the codon option is altered. 

## Using Snakemake

VaxPress can be run in a GCP (Google Compute Platform) VM instance using Snakemake. The `Snakefile` is provided in the `snakemake` directory. The `gcloud-setup.sh` file contains the setup process for the google compute platform. The snakemake depends on several configurations. 

Inside a directory which looks like this:
```
  test
   ├──  sequence.fasta
   └──  config.json
```
do

```bash
cd test
snakemake -p --snakefile /path/to/Snakefile \ 
    --config \
     input=sequence.fasta \
     vaxpress-config=parameters.json\
     output-dir='wow_it_works'\
     instance-name='vaxpress-test'
```

the results will go inside the `wow_it_works` directory inside the `test` directory.

# Citing VaxPress

If you employed our software in your research, please
kindly reference our publication:

Ju, Ku, and Chang (2023) Title. Journal. Volume. (in preparation)

# License

VaxPress is distributed under the terms of the [MIT License](LICENSE.txt).

LinearFold and LinearDesign are licensed for non-commercial use
only. If considering commercial use, be cautious about using options
such as `--lineardesign` and `--folding-engine linearfold`.
