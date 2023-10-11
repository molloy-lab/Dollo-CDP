Example for Myotis
------------------

We now show how to run Dollo-CDP using a retroelement insertion presence/absence data set from [Korstian et al. (2022)](https://doi.org/10.3390/genes13030399). 

First, we downloaded the data set from [this Github repository](https://github.com/jkorstia/retrophylogenomic_tools/blob/main/myotis_ves.nex). Second, we updated the Nexus format because our current codebase is somewhat picky. An example of the format, with the first 50 characters, is below. Note that missing values should be indicated with a `?`.

```
#NEXUS

Begin data;
	Dimensions ntax=11 nchar=50;
	Format datatype=standard gap=-;
	Matrix
Aust 11111011110111110111011101111111111011101111111101
Bran 00001011000001011110101000100000101100010010000000
Cili 11111111011011110011110110111101100010110111111101
Davi 11010000000001100000100101001101001110101010000000
Luci 11111111111111111111111111111111111111111111111111
Occu 11111101011101110110111111111101111011111110110110
Sept 01110101110111111111111111101111011111111111111111
Thys 11000010100110011110001110100100101001001001001011
Veli 11110111111001111111111011111111011111111111111110
Vive 10011110101100100000101101011011000110111010100011
Yuma 11111111101111101101010111011111001111111101010101
	;
End;
```

To run Dollo-CDP on the Myotis data set, use the following command:

```
../src/dollo-cdp \
    -i myotis_ves.nex \
    -g Davi \
    -o dollo-cdp-myotis-ves.tre
```

The output tree should be
```
((Davi),((Bran),(((Sept),((Cili),(Thys))),(((Vive),((Luci),(Occu))),((Aust),((Yuma),(Veli)))))));
```
This tree has a Dollo score of `11618` and is topologically equivalent to the Dollo tree presented in Figure 2A by Korstian et al.

Because this data set has only 11 taxa, we can also analyze this data set using branch-and-bound. We performed this analysis using [PAUP*](https://paup.phylosolutions.com), although you might also try the [Dollop function distributed with Phylip](https://evolution.genetics.washington.edu/phylip/doc/dollop.html).

After downloading PAUP* to this directory, use the following commands, making sure to use the correct binary name:
```
gunzip paup4a168_osx
chmod a+x paup4a168_osx
./paup4a168_osx -n paup-dollo-bnb-myotis-ves.nex 
```
The resulting tree is topologically equivalent to the one recovered by Dollo-CDP. Note that PAUP* reports a Dollo score of `20939`, which differs from our score because we only count losses (`1 -> 0`) and do not count gains (`0 -> 1`). You can score the PAUP tree using the command:
```
../src/dollo-cdp \
    -q paup-dollo-bnb-myotis-ves-strict.tre \
    -i myotis_ves.nex 
```
The output score is `11618`.

To analyze data sets with larger numbers of taxa, you may wish to speed-up Dollo-CDP by further constraining the space of allowed solutions. This can be achieved by performing heuristic search with PAUP* and then giving the resulting trees to Dollo-CDP as constraints. Let's create a constraints by building a starting tree via random taxon addition and then performing TBR moves, saving the 100 best scoring trees.
```
./paup4a168_osx -n paup-dollo-hsearch-myotis-ves.nex 
```
Now let's give these trees as constraints to Dollo-CDP.
```
../src/dollo-cdp \
    -i myotis_ves.nex \
    -g Davi \
    -t paup-dollo-hsearch-myotis-ves-all.trees \
    -o dollo-cdp-plus-hsearch-constraints-myotis-ves.tre
```
This gives us the same tree as the other two analyses with Dollo parsimony.

A copy of the PAUP* user manual is available [here](https://phylosolutions.com/paup-documentation/paupmanual.pdf); also see [https://rothlab.ucdavis.edu/genhelp/paupsearch.html](https://rothlab.ucdavis.edu/genhelp/paupsearch.html).
