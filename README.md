# PlasEval

PlasEval is a tool aimed at evaluating the accuracy and at comparing methods for the problem of **plasmid binning**.

## Plasmid binning

Plasmid binning aims at detecting, from the draft assembly of a single isolate bacterial pathogen, groups of contigs assumed eachto originate from a plasmid present in the sequenced isolate. Recent methods for plasmid binning include [PlasBin-Flow](https://github.com/cchauve/PlasBin-flow), [MOB-recon](https://github.com/phac-nml/mob-suite) and [gplas](https://gitlab.com/sirarredondo/gplas).

The result of applying a plasmid binning method to a draft assembly is a set of unordered groups of contigs, called **plasmid bins**.

## Plasmid binning evaluation and comparison

### Evaluation

In the **evaluation** mode of PlasEval, given a set of *predicted plasmid bins* resulting from a plasmid binning tools and a *ground truth* set of plasmid bins, where each true plasmid bin contains all the contigs that belong to one of the plasmid present in the sequenced isolate, PlasEval computes three statistics, the *precision*, the *recall* and the *F1-score* of the predicted plasmid bins with respect to the ground truth.

For a group $`X`$ of contigs, we denote by $`L(X)`$ the cumulated length of the contigs in $`X`$.
For a predicted plasmid bin $`P`$ and a ground truth plasmid bin $`T`$, we define the *overlap* between $`P`$ and $`T`$, denoted by $`o(P,T)`$, as the cumulated length of the contigs presents in both $`P`$ and $`T`$.
Given a set $`A`$ of predicted plasmid bins and a set $`B`$ of ground truth plasmid bins, we define the precision $`p(A,B)`$ and recall $`r(A,B)`$ as follows:

```math
p(A,B) = \frac{\sum\limits_{P\in A} \max\limits_{T\in B} o(P,T)}{\sum\limits_{P\in A}L(P)}, r(A,B) = \frac{\sum\limits_{T\in B} \max\limits_{P\in A} o(P,T)}{\sum\limits_{T\in B}L(T)}.
```

The F1-score is the arithmetic mean of the precision and recall

```math
F_1(A,B) = 2\frac{{p}(A,B){r}(A,B)}{{p}(A,B)+{r}(A,B)}.
```

### Comparison

In the **comparison** mode of PlasEval, given two sets of *predicted plasmid bins* (either resulting from two plasmid binning tools or from a plasmid binning tool and a ground truth), PlasEval computes three statistics, PlasEval computes a *dissimilarity measure* that indictaes how much both sets of predicted plasmid bins are in agreement.
The full details of the dissimilarity measure computed by PlasEval are available in the preprint *[PlasEval: a framework for comparing and evaluating plasmid binning tools](https://link.springer.com/article/10.1186/s12859-024-05941-0)* and we describe below its general principle.
The dissimilariy value is the sum of 4 terms accounting respectively for

- *extra contigs*: the contigs present in $`A`$ but not in $`B`$;
- *missing contigs*: the contigs present in $`B`$ but not in $`A`$;
- *splits*: splitting of the plasmid bins of $`A`$ to obtain a set of intermediate bins defined as the intersection of $`A`$ and $`B`$;
- *joins*: joining the splitted plasmid bins into the plasmid bins of $`B`$.

Each of these components incur a *cost*, parameterized by a parameter $`\alpha \in [0,1]`$:

- each extra or missing contig $`c`$ results in a cost $`\ell(c)^\alpha`$, where $`\ell(c)`$ is the length of $`c`$;
- each split of a plasmid bin $`P`$ into two smaller bins $`P',P''`$ results in a cost $`\min(L(P')^\alpha,L(P'')^\alpha)`$;
- each join of two intermediate plasmid bins $`P',P''`$ into  larger bin $`P`$ results in a cost $`\min(L(P')^\alpha,L(P'')^\alpha)`$.

So with $`\alpha=0`$, the length of contigs is not accounted for and only set-theoretic operations define the dissimilarity value, while with $`\alpha=1`$ it is fully accounted for.
By default $`\alpha=0.5`$.

In the case where some contigs are *repeated*, i.e. a contig appears in more than one plasmid bin of $`A`$ and/or $`B`$, a *branch-and-bound* algorithm computes the pairing between repeats contigs in $`A`$ and in $`B`$ that results in the minimum dissimilarity value.

Finally the dissimilarity obtained as described above is *normalized* into a value in $`[0,1]`$ by dividing it by

```math
\sum\limits_{P\in A}\sum\limits_{c\in P} \ell(c)^{\alpha} + \sum\limits_{Q\in B}\sum\limits_{c\in Q} \ell(c)^{\alpha}.
```

## Installation

```sh
pip install -r requirements.txt
```

## Usage

### Input: plasmid bins file

The main input in both modes of PlasEval are the two plasmid bins files.
A plasmid bins file is a TSV file that describes a set of plasmids bins in a forma where each row contains three pieces of information: the identifier of a plasmid bin, the identifier of a contig that belongs to this plasmid bin and the length of the contig.
The file should have a header row with column names `plasmid`, `contig`, `contig_len`.
An example is provided below that describes the set of plasmid bins where bin `P1` contains contigs `C1,C2` and bin `P2` contains contigs `C1`, `C3` and `C4`.

```txt
plasmid contig  contig_len
P1 C1  2000
P2 C3  3000
P1 C2  2000
P2 C1 2000
P2 C4 2000
```

If a contig appears in several copies in a plasmid bin, the evaluation mode only accounts for one copy of the contig. The comparison mode can account for multiple copies of a contig in the same bin.

### Input: numeric parameters

In both evaluation and comparison mode, PlasEval takes an extra optional parameter `min_len`: every contig of length below the value `min_len` is discarded from both sets of considered plasmid bins. This parameter is useful in comparison mode in the case of plasmid bins sets that contain many short repeated contigs, which can result in the branch-and-bound algorithm taking a long time to complete.

The comparison mode uses two more parameters. Firstly, the value of $\alpha$ can be passed as a parameter `p`, although by default it takes value $0.5$. Secondly, the maximum number of recursive calls used in the branch and bound can also be set by the user. If the number of recursive calls is exceeded, the comparison is stopped. In such instances, the comparison mode can be rerun with a higher length threshold. The default value for the maximum number of recursive calls (`max_calls`) is $10000000$.

### Runing `eval` or `comp` modes

1. The following command is used for the evaluation mode:

   ```sh
   python plaseval.py eval --pred PREDICTED_BINS_TSV --gt GROUNDTRUTH_BINS_TSV --out_file OUT_FILE --log_file LOG_FILE (--min_len LEN_THRESHOLD)
   ```

   Where `pred` and `gt` are TSV files, with the set of predicted and ground truth plasmid bins respecitvely. `out_file` is the path to the output file. The integer length threshold `min_len` can be provided as an optional parameter.

2. In addition to the two plasmid bin files, the comparison mode also takes the path to the log file as input.
The following command is used for the comparison mode:

   ```sh
   python plaseval.py comp --l LEFT_BINS_TSV --r RIGHT_BINS_TSV --out_file OUT_FILE --log_file LOG_FILE (--min_len LEN_THRESHOLD --p ALPHA --max_calls MAX_RECURSIVE_CALLS)
   ```

   Where `LEFT_BINS_TSV` and `RIGHT_BINS_TSV` are TSV files, each with one set of plasmid bins. `out_file` is the path to the output file while `log_file` is the path to the log file. The parameters `min_len`, `p` and `max_calls` are optional.

### Output

1. The output file is TSV file with the following columns:
   1. `Level`: Precision and recall statistics are computed for individual bins (`Individual` level). The average precision, recall and F1 statistics are then computed for the overall sample (`Overall` level).
   2. `Statistic`: This columns list the type of statistic (`Precision`, `Recall` or `F1`).
   3. `Bin`: For the `Individual` level, we compute the precision for each predicted plasmid bin and the recall for each ground truth plasmid bin. The identity of the bin (predicted or ground truth) is given in this column. Note that this column is empty for `Overall` sample statistics.
   4. `Unwtd_Stat`: Contig-level statistics for an individual bin or for the overall sample.
   5. `Wtd_Stat`: Basepair-level statistics for an individual bin or for the overall sample.
   6. `Unwtd_Match`: The best bin from the opposite side matched to the bin in question (from the 'Bin' column) according to contig-level statistics. For `Precision`, this column will have the ground truth bin that best matches the predicted bin from the 'Bin' column. For `Recall`, this column will have the predicted bin that best matches the ground truth bin from the 'Bin' column. Note that this column is empty for `Overall` sample statistics.
   7. `Wtd_Stat`: The best bin from the opposite side matched to the bin in question (from the 'Bin' column) according to basepair-level statistics. For `Precision`, this column will have the ground truth bin that best matches the predicted bin from the 'Bin' column. For `Recall`, this column will have the predicted bin that best matches the ground truth bin from the 'Bin' column. Note that this column is empty for `Overall` sample statistics.

2. The output file for the compare mode contains the following information:
   1. `Total_ctg_length`: Cumulative length of contigs present in at least one of set of plasmid bins.
   2. `Total_ctg_length_alpha`: Cumulative dissimilarity cost of all contigs from (a).
   3. `Cuts`: Cost of cuts: cost of splitting bins from the first set of plasmid bins.
   4. `Joins`: Cost of joins: cost of splitting bins from the second set of plasmid bins.
   5. `Extra_ctgs`: Cumulative length of contigs present only in the first set.
   6. `Missing_ctgs`: Cumulative length of contigs present only in the second set.
   7. `Dissimilarity`: Dissimilarity score

The compare mode also provides a log file with some other details related to the comparison algorithm. These include the maximum number of matchings possible, the time taken to execute the method, the number of recursive function calls made during the comparison and finally the actual matching between contigs of both sets of plasmid bins that yields the dissimilarity score in the output file described above.

### Examples

A few toy examples to demonstrate the use of PlasEval have been provided in the examples directory.

1. The following command evaluates predicted plasmid bins (`pred_bins_1.tsv`) against the ground truth plasmid bins (`gt_bins_1.tsv`). Details of the evaluation output will be printed to the file `P1G1_eval.out`.

   ```sh
   cd src/
   python plaseval.py eval --pred ../examples/input/pred_bins_1.tsv --gt ../examples/input/gt_bins_1.tsv --out_file ../examples/output/P1G1_eval.out --log_file ../examples/output/P1G1_eval.log
   ```

2. The following command evaluates predicted plasmid bins (`pred_bins_1.tsv`) against the other set of ground truth plasmid bins (`gt_bins_2.tsv`), only considering contigs above $1000$ bp. Details of the evaluation output will be printed to the file `P1G2_eval.out`.

   ```sh
   cd src/
   python plaseval.py eval --pred ../examples/input/pred_bins_1.tsv --gt ../examples/input/gt_bins_2.tsv --min_len 1000 --out_file ../examples/output/P1G2_eval.out
   ```

3. The following command will compare the sets of plasmid bins in `pred_bins_1.tsv` and `gt_bins_1.tsv`. The details of the dissimilarity between the two sets will be output to `P1G1_comp.out` while the log file `P1G1_comp.log` contains the miscellaneous detials.

   ```sh
   cd src/
   python plaseval.py comp --l ../examples/input/pred_bins_1.tsv --r ../examples/input/gt_bins_1.tsv --out_file ../examples/output/P1G1_comp.out --log_file ../examples/output/P1G1_comp.log
   ```

4. The following command will compare the sets of plasmid bins in `pred_bins_1.tsv` and `pred_bins_2.tsv`, neither of which contains the ground truth. The cost parameter $\alpha$ is set at $0.75$. The details of the dissimilarity between the two sets will be output to `P1P2_comp.out` while the log file `P1P2_comp.log` contains the miscellaneous detials for the comparison.

   ```sh
   cd src/
   python plaseval.py comp --l ../examples/input/pred_bins_1.tsv --r ../examples/input/pred_bins_2.tsv --p 0.25 --out_file ../examples/output/P1P2_0.75_comp.out --log_file ../examples/output/P1P2_0.75_comp.log
   ```
