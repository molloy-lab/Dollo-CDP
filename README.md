Dollo-CDP
=========
Dollo-CDP is a method for estimating phylogenies from character data under the Dollo parsimony criterion score. 
The input characters must be ordered, with 0 representing the ancestral state and 1 representing the mutated or derived state. 
Under the Dollo assumption, there can be exactly one forward mutation (0->1) per character, however, there can be any number of backward mutations (1->0); this assumption is popular for analyzing low-homoplasy characters evolving within a population genetics framework (i.e., where there is ILS) as well as in tumor phylogenetics.
Our approach is novel in that it is guaranteed to return an optimal solution to the Dollo parsimony problem that obeys the set of constraints (clades) given as input.
These constraints can be generated automatically from the character data or via a previous heuristic search, as described in [this example](example/README.md).

To learn more about Dollo-CDP, check out the related paper in [Algorithms Mol Biol  / WABI 2023](https://doi.org/10.1186/s13015-023-00249-9).

Usage
-----

To build, Dollo-CDP use commands:
```
git clone https://github.com/molloy-lab/Dollo-CDP.git
cd Dollo-CDP/src
make
```
Note: On Linux, we successfully compiled with gcc version 8.5.0 and version 9.3.0. On Mac OS X, we successfully compiled with Apple clang version 14.0.3; we were unable to compile with gcc installed via homebrew, unfortunately. The former requires Apple command line tools to be installed. This can be done with the following commands
```
# sudo rm -rf /Library/Developer/CommandLineTools
xcode-select --install
```
and then following the pop-up.

Alternatively, you could download binaries in a release. In either case, before running Dollo-CDP, you must download ASTRAL and extract the zip folder into the src directory:
```
git clone https://github.com/smirarab/ASTRAL.git
mv ASTRAL tmp-ASTRAL
unzip tmp-ASTRAL/Astral.*.zip
rm -rf tmp-ASTRAL
```

To run Dollo-CDP, use command:
```
./dollo-cdp -i <input nexus file> -o <output file name>
```
We recommend working through [this example](example/README.md).

We also recommend checking out the Dollo-CDP usage options with this command:
```
./dollo-cdp -h
```
The output should be
```
Dollo-CDP version 1.0.0
COMMAND: ./dollo-cdp 
===================================== Dollo-CDP =====================================
Dollo-CDP is a program that solves the large Dollo parsimony problem for binary
characters (missing values allowed) within a clade-constrained version of tree space.

USAGE for Large Dollo problem:
./dollo-cdp -i <input characters file> -g <outgroup name> -o <output file>

USAGE for small Dollo parsimony problem:
./dollo-cdp -i <input characters file> -q <input species tree>

OPTIONS:
[-h|--help]
        Prints this help message.
(-i|--input) <input characters file>
        Name of file containing input characters in nexus format (required)
(-x|-g|--outgroup) <outgroup name>
        Comma separated list of outgroup taxa used to root solution space
[(-t|--trees) <input trees file>]
        Name of file containing trees in newick format for constructing solution
        space with ASTRAL. If no file is specified, characters will be treated as
        rooted trees with at most one internal branch.
[(-q) <input species file>]
        Name of file containing species trees in newick format
[(-k)]
        Write characters as rooted trees with one internal branch and exit
[(-o|--output) <output file>]
        Name of file for writing output species tree (default: stdout)

Contact: Post issue to Github (https://github.com/molloy-lab/Dollo-CDP)
        or email Junyan Dai (jdai1234@umd.edu) & Erin Molloy (ekmolloy@umd.edu)

If you use Dollo-CDP in your work, please cite:
  Dai, Rubel, Han, Molloy, 2024, "Dollo-CDP: a polynomial-time algorithm for
  the clade-constrained large Dollo parsimony problem", Algorithms for Molecular
  Biology, https://doi.org/10.1186/s13015-023-00249-9.
====================================================================================
```

