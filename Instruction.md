# Instructions
* [Datasets](#datasets)
  * [GIFGIF](#gifgif)
  * [CHESS](#chess)
  * [NASCAR](#nascar)
* [Numerical Experiments](#experiments)
  * [Bradley-Terry model - GIFGIF dataset](#bt)
    * [Full Dataset](#btfulldataset)
    * [Sample Dataset](#btsampledataset)
  * [Rao-Kupper model - Chess dataset](#rk)
    * [Full Dataset](#rkfulldataset)
    * [Sample Dataset](#rksampledataset)
  * [Plackett-Luce model - NASCAR dataset](#pl)

## Datasets

### [GIFGIF](http://lucas.maystre.ch/gifgif-data)
This dataset contains user evaluations of digital images by paired comparisons with respect to different metrics, such as amusement, content, and happiness. The dataset was collected through an online web service by the MIT Media Lab as part of the PlacePulse [project](http://gifgif.media.mit.edu/). This service presents the user with a pair of images and asks to select one that better expresses a given metric, or select neither. The dataset contains 1,048,576 observations and covers 17 metrics. We used this dataset to evaluate convergence of MM algorithms for **Bradley-Terry** model of paired comparisons. We did this for the three aforementioned metrics.

### [Chess](https://www.kaggle.com/c/chess/data)
This dataset contains game-by-game results for 65,030 matches among 8,631 chess players. The dataset was used in a Kaggle chess ratings [competition](https://www.kaggle.com/c/chess). Each observation contains information for a match between two players including unique identifiers of the two players, information about which one of the two players played with white figures, and the result of the match, which is either win, loss, or draw. This dataset has a large degree of sparsity. We used this dataset to evaluate convergence of the **Rao-Kupper** model of paired comparisons with ties.

### [Nascar](http://personal.psu.edu/drh20/code/btmatlab/)
This dataset contains auto racing competition results. The source of this dataset is [here](http://personal.psu.edu/drh20/code/btmatlab/). Each observation is for an auto race and contains the ranking of drivers in increasing order of their finish times in the race. We used this dataset to evaluate convergence times of MM algorithms for the **Plackett-Luce** ranking model.

## Numerical Experiments <a name="experiments"></a>
### Bradley-Terry model - GIFGIF dataset <a name="bt"></a>

* [GIFGIF](http://lucas.maystre.ch/gifgif-data) dataset contains outcomes for 17 different metrics (happiness, pride, fear, amusement, contempt, ...)
* In our experiments, we used three metrics (amusement, contempt, and happines)
* The experiments are conducted by using full dataset for the three metrics, and using a sample from the strongly connected coponent.

#### Step-by-step instructions:

##### Full dataset <a name="btfulldataset"></a>

**Step 1**. Use 'read_GIFGIG_data.R' to create the winning matrix ('w') and the comparison matrix ('A'); this also 
computes 'd(M)' (m.max) and 'a(M)' (second smallest eigenvalue).

**Step 2**. In 'MMandGDbt.R', we can compute the number of iterations until eps-convergence for both GD and MM algorithm. 
For GD algorithm, we use function

```r
Bayesian.GD(init,tiny,alpha,beta,w,m.max)
```

and for MM algorithm, we use function

```r
Bayesian.MM(init,tiny,alpha,beta,w)
```
**Step 3**. For both functions, 'init' is the initial vector value of parameters, 'tiny' is eps. 'w' is the winning matrix and 'm.max' is the maximum number of paired comparison that we computed from 'read_GIFGIF_data.R'. Each of these function will return a matrix called ('Iterationmatrix'). Iterationmatrix(i,j) is the updated parameter value for j-th node at the i-th iteration. The total number of iterations until eps-convergence can be obtained for both GD and MM case using the following code

```r
#For GD case
dim(Bayesian.GD(init,tiny,alpha,beta,w,m.max)$Iterationmatrix)[1]

#For MM case
dim(Bayesian.MM(init,tiny,alpha,beta,w)$Iterationmatrix)[1]
```
**Step 4**. Repeat the experiment for beta=0.01, 0.1, 1, 10 respectively and then record the number of iterations for each pair of alpha and beta.

##### Sample dataset <a name="btsampledataset"></a>

**Step 1**. Use 'construct_strong_connected_graph.R' to extract the strongly connected sample graph from full dataset for each metric. 

**Step 2**. In 'construct_strong_connected_graph.R', we can compute the winning matrix ('w0') and the comparison matrix ('AA'); this also 
computes 'd(M)' (m.max) and 'a(M)' (second smallest eigenvalue) for the sample graph.

**Step 3**. Follow **step 2, 3 and 4** from [Full dataset](#btfulldataset) to compute the number of interations for each pair of alpha and beta.


## Rao-Kupper model - Chess dataset <a name="rk"></a>

* [Chess](https://www.kaggle.com/c/chess/data) contains results for 65,350 matches among 8631 chess players.

#### Step-by-step instructions:

##### Full dataset <a name="rkfulldataset"></a>

**Step 1**. Use 'read_Chessrating_data.R' to create the winning matrix ('s') and the comparison matrix ('A'); this also 
computes 'd(M)' (m.max) and 'a(M)' (second smallest eigenvalue).

**Step 2**. In 'MMandGD_Rao_Kupper.R', we can compute the number of iterations until eps-convergence for both GD and MM algorithm. 
For GD algorithm, we use function

```r
Bayesian.GD(init,tiny,alpha,beta,s,m.max,c)
```

and for MM algorithm, we use function

```r
Bayesian.MM(init,tiny,alpha,beta,s,c)
```
**Step 3**. For both functions, 'init' is the initial vector value of parameters, 'tiny' is eps and 'c'= sqrt(2). 's' is the winning matrix and 'm.max' is the maximum number of paired comparison that we computed from 'read_Chessrating_data.R'. Each of these function will return a matrix called ('Iterationmatrix'). Iterationmatrix(i,j) is the updated parameter value for j-th node at the i-th iteration. The total number of iterations until eps-convergence can be obtained for both GD and MM case using the following code

```r
#For GD case
dim(Bayesian.GD(init,tiny,alpha,beta,s,m.max,c)$Iterationmatrix)[1]

#For MM case
dim(Bayesian.MM(init,tiny,alpha,beta,s,c)$Iterationmatrix)[1]
```
**Step 4**. Repeat the experiment for beta=0.01, 0.1, 1, 10 respectively and then record the number of iterations for each pair of alpha and beta.

##### Sample dataset <a name="rksampledataset"></a>

**Step 1**. Use 'construct_strong_connected_graph.R' to extract the strongly connected sample graph from full dataset. 

**Step 2**. In 'construct_strong_connected_graph.R', we can compute the winning matrix ('s0') and the comparison matrix ('AA'); this also 
computes 'd(M)' (m.max) and 'a(M)' (second smallest eigenvalue) for the sample graph.

**Step 3**. Follow **step 2, 3 and 4** from [Full dataset](#rkfulldataset) to compute the number of interations for each pair of alpha and beta.


## Plackett-Luce model - NASCAR dataset <a name="pl"></a>

* [Nascar](http://personal.psu.edu/drh20/code/btmatlab/) contains auto racing competition results for 36 races with 83 players in total.

#### Step-by-step instructions:

**Step 1**. Use 'read_NASCAR_data.R' to create the  matrix ('playersmatrix') with playersmatrix(i,j) represents as the player's ID for i-th race and j-th ranked player; this also computes the comparison matrix ('A'), 'd(M)' (m.max) and 'a(M)' (second smallest eigenvalue).

**Step 2**. In 'MMandGD_PL.R', we can compute the number of iterations until eps-convergence for both GD and MM algorithm. 
For GD algorithm, we use function

```r
Bayesian.GD(init,tiny,alpha,beta,playersmatrix,n)
```

and for MM algorithm, we use function

```r
Bayesian.MM(init,tiny,alpha,beta,playersmatrix,n)
```
**Step 3**. For both functions, 'init' is the initial vector value of parameters, 'tiny' is eps. 'playersmatrix' is the matrix ('playersmatrix') that we computed from 'read_NASCAR_data.R'. Each of these function will return a matrix called ('Iterationmatrix'). Iterationmatrix(i,j) is the updated parameter value for j-th node at the i-th iteration. The total number of iterations until eps-convergence can be obtained for both GD and MM case using the following code

```r
#For GD case
dim(Bayesian.GD(init,tiny,alpha,beta,playersmatrix,n)$Iterationmatrix)[1]

#For MM case
dim(Bayesian.MM(init,tiny,alpha,beta,playersmatrix,n)$Iterationmatrix)[1]
```
**Step 4**. Repeat the experiment for beta=0.01, 0.1, 1, 10 respectively and then record the number of iterations for each pair of alpha and beta.
