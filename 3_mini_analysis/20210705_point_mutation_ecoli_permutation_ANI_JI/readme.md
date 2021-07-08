### Analysis plan for JI vs ANI in E.coli permutation model

---

#### Aim: 

Correlates average nucleotide identity (ANI) and Jaccard index (JI) by k-mers in range of k values.



#### Parameter:

1. Input genome (default: e.coli for now), the baseline
2. k range (default: 15-60 by 5), length for k-mer
3. p range (default: 0.05-0.9 by 0.05), simple point mutation
4. simulation size m (default: 1000)



#### More details to confirm:

1. independent & simple point mutation: 
   1. generate n numbers ~ Unif(0,1) and use p as the cutoff for 0/1
   2. mutation target: equally to the rest 3 letters
2. Expectation of JI/CI: that's empirical mean or median?



#### Analysis plan:

1. data preparation:
   1. transfer input genome -> a long string S
2. loop through mutation rate p in p range, for each p value:
   1. generate m mutated strings
   2. for every mutated string S', get the following between S and S':
      1. ANI: close to p, for QC purpose only in case of small n
      2. JI for all k values
      3. CI for all k values
   3. plot the distribution of JI, CI by k value. Get the **mean/median** to compare with the theoretical expectation.



#### Output:

1. many line/box plots of ANI vs JI: x axis is k, y axis is JI.
2. many tables records the expectation of JI/CI 
