# Dollo-CDP

This repository corresponds to the method described in this paper: ``Dollo-CDP: Leveraging constraints plus dynamic programming for the large Dollo parsimony problem.''

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
COMMAND: ./dollo-cdp -h 
Not enough arguments given!

=================================== Dollo-CDP ===================================
This is version 1.0.0 of Dollo-CDP, a program that solves the large Dollo
parsimony problem within a clade-constrained version of the solution space.

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
================================================================================
```

