# Description of test cases

The folder `test_cases/input` contains several "test" plasmid bins. Pairs of these test bins have been to compared against each other to test different functionalities of PlasEval. Given below is the list of comparisons made along with the significance of each comparison.

1. Bins 0 vs 1: Test bins 0 is empty. This tests if an empty bin is appropriately dealt with. The dissimilarity score should be 1.

2. Bins 0 vs 0: Test bins 0 is empty. This tests the case where both sets of bins are empty. The dissimilarity score should be 0.

3. Bins 1 vs 1: Tests two exactly same bins. The dissimilarity score should be 0.

4. Bins 1 vs 1 with alpha = 0: Tests two exactly same bins but with alpha = 0. The dissimilarity score should still be 0.

5. Bins 2 vs 2: Tests two exactly same bins with copies. The dissimilarity score should be 0.

6. Bins 2 vs 2 with alpha = 0: Tests two exactly same bins with copies but with alpha = 0. The dissimilarity score should still be 0.

7. Bins 1 vs 2: Tests two bins with same contig families but extra copies (in test bins 2).

8. Bins 3 vs 4: Tests two bins with same contig families but extra copies on both sides.

9. Bins 1 vs 5: Tests two bins with same contig families but extra copies in the same bin (in test bins 5).

10. Bins 1 vs 6: Completely different bins. The dissimilarity score should be 1.

11. Bins 1 vs 7: Simple split. Only splits cost should be incurred.

12. Bins 1 vs 8: Simple join. Only joins cost should be incurred.

13. Bins 1 vs 9: Splits + missing contigs.

14. Bins 1 vs 10: Splits + extra contigs.

15. Bins 1 vs 11: Splits + missing copies.

16. Bins 7 vs 8: Split into multiple bins.

17. Bins 8 vs 7: Join from multiple bins.

18. Bins 7 vs 12: Multiple splits and joins.

19. Bins 9 vs 14: Multiple splits and joins with extra contigs on both sides.

20. Bins 2 vs 11: Extra copies on both sides. Illustrates branch-n-bound functionality.

21. Bins 2 vs 13: Extra copies on both sides. Illustrates branch-n-bound functionality.

NOTE: Tests 20 and 21 were designed to test cases with splits and joins with extra copies on both sides. However, they ended up without any splits (both tests 20, 21) or joins (test 20). Interestingly, the algorithm matched the contig copies in such a way that the splits / joins were not required to transform one set of bins into the other.
