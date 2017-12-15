---
title: "A Memo on Statistical Computing"
author: "Jiayu Wu"
output: 
       html_document:
              toc: true
              theme: united
---

# Preface
This is a memo on the course Statistical Computing (202A, 2017Fall, UCLA). It is mainly about the intuition and implementation of fundamental statistical learning methods, using R as primary coding lanuage facilitated by other computing languages and tools. 

The content is divided into three parts. The first part reviews and compares different alorithms in terms of their intuition, computation tricks and efficiency. The second part summarizes experiences with various computing tools in this course. The last part is the coding implementation with comments. To view codes of a particular algorithm in part one, click [*View Code*] in text.

# Algorithm

## OLS Regression

We start from the classical linear regression model from the perspective of statistical learning. The process of estimating the unknown model parameter $beta$ using known dataset is learning from training data, for the purpose of explanation and/or prediction. With least square estimation minimizing RSS we can derive an estimator for $beta$, yet the computation could be complicated.

We write compuiter programs to do the complex computation for us, and we are looking for more efficient algorithm so that the computation is still fast with a large quantity of data.

### Gauss Jordan Elimination [[*View Code*]](#GJ)
Gauss Jordan Elimination is a method to find inverse of a matrix based on row operations, and two calculus tricks make it more powerful for linear regression:   
1) Partitioned the matrix to derive a matrix version of Gauss-Jordan;   
2) Contruct a cross product matrix $[X \; Y]^T[X \; Y]= \begin{bmatrix} X^{T}X & X^{T}Y \\ Y^{T}X & Y^{T}Y \end{bmatrix}$     

Thus we get
$\begin{bmatrix} X^TX & X^{T}Y & | & I_1 & 0 \\ Y^{T}X & Y^{T}Y & | & 0 & I_2 \end{bmatrix}\stackrel{GJ[1:m]}{\to}\begin{bmatrix} I_1 & \hat{\beta} & | & Var(\hat{\beta})/\sigma^2 & 0 \\ 0& RSS & | & -\hat{\beta}^T & I_2 \end{bmatrix}$.   
which only requires m (number of columns of design matrix X) times of simple calculation to get the estimation and the standard error.

### Sweep Operator [[*View Code*]](#swp)
Sweep Operator is a more compact version of Gauss-Jordan. It uses the same mean of calculation but save space by omitting identity matrix in Gauss Jordan:

$\begin{bmatrix} X^TX & X^{T}Y \\ Y^{T}X & Y^{T}Y & \end{bmatrix}\stackrel{SWP[1:m]}{\to}\begin{bmatrix} -Var(\hat{\beta})/\sigma^2 & \hat{\beta} &  \\ \hat{\beta}^T & RSS \end{bmatrix}$.   

### QR Decomposition [[*View Code*]](#qr)
QR Decomposition is an orthogonal rotation of the data in a high-dimensional space. This rotation can be understanded as a change of perpective, through which we view the data from a different basis. In this basis X becomes an upper triangular matrix R whose inverse is easy to compute. 
To perfrom this transformation, we apply Householder reflection. The key insight is that orthogonal transformation preserves the length of vectors, i.e. the norm of each column of R equals that of X. Also we know that R is upper triangular, hence the transformation can be identified. Note that the sign of element in R should be chosen as the opposite of corresponding element in X for numerical stability.

### Principal Component Analysis [[*View Code*]](#pca)
QR method can also be used to perform principal component analysis for dimension reduction. The idea of PCA is to rotate the high-dimensional space spanned by X (note that X should be centralized first) such that the variation of the observations is small on some of the axises from the new basis, and these axis can be omitted without a big loss of information from the data, thus we achieve reduction of dimension. 

PCA can be used in both supervised and unsupervised learning. Regression is a typical example of the former with a dependent variable of interest (y), where PCA can be used to solve the multicollinearity problem. In unsupervised learning, PCA is also powerful in unsupervised learning to explore the structure of X by condense the information to the most representative features. 

Realization of PCA involves eigen decomposition of the variance-covariance-matrix onto a orthogonal basis. Such a symmetrix matrix can be diagnonalized as $\Sigma=Q\Lambda Q^T$ where Q is an orthogonal matrix (can be computed with QR method) and $\Lambda$ is a diagonal matrix, it can be interpretated as that $\Sigma$ is $\Lambda$ in base Q. To find the Q and $\Lambda$ in question, we start from any square matrix V, then orthogonalize it and multiply it by $\Sigma$ to update the matrix V, by iteration it will converge to the product of Q and $Lambda$. 

## Regularized Learning

Regularized learning is a way to avoid overfitting, that is to avoid interpreting noise $(X^TX)^{-1}X^T\epsilon$. In classical modeling, we test overfitiing with hypothesis testing, and exclude $\beta_k$ falls out of $2\sigma$, it is refereed to as hard-thresholding.
From the perspective of statistical learning, complex model with more parameters always has smaller training error, but it would increase testing error. The idea of regularization is to add bias to the estimators $\beta_k$ towards 0, so we tend to choose a simpler model with less parameters. Then it is critical to trade off between bias and variance, i.e., to find a ideal $\lambda$. We use cross-validation for this purpose: split the data for training and testing, then find the point minimize their difference.

### Ridge Regression [[*View Code*]](#ridge)
Ridge Regrssion minimize the following loss function:

$\hat{\beta}_{\lambda}  = \arg\min_{\beta} \left[ \|Y - X \beta\|_{\ell_2}^{2}  + \lambda \|\beta\|_{\ell_2}^2\right] = (X^{T}X + \lambda I_p)^{-1}X^{T}Y$

It can be considered as add to linear regression a penalty term on $\beta$ to prefer smaller estimator. When $(X^TX)^{-1}=I$, we have $\hat{\beta}_{k,ridge}=\frac{\hat{\beta}_{k,ls}}{1+\lambda}$. So ridge can be viewed as a shrinkage towards 0.

In computation implementation, we could build on least square methods with the idea of "suedo data". That is, we assume the design matrix to be $(X \; \Lambda=diag[\sqrt{\lambda}])^T$, response to be $(Y\;0)^T$. Intercept should not be pernalized because it does not have an interpretation of correlation.

### Spline Regression [[*View Code*]](#spline)
Spline Regression is an example of use of ridge regression. The idea is to use piecewise linear regression to approximate non-linear function. The more the knots we set for the piecewise regression, the more the hidden variables and parameters. When the number of knots (number of parameters) is close or greater than the number of observations, the model could even interpret each obervation one by one, which is the overfitting problem.

By add a regularization term, we bias $\hat{\beta_k}$ towards 0 in order to get a smooth function. This example shows that, ridge regression as shrinkage is especially useful for situation where p>n.

### Lassco Regression [[*View Code*]](#lassco)
Lassco Regression minimize the following loss function:

$\hat{\beta}_{\lambda}  = \arg\min_{\beta} \left[  \frac{1}{2} \|Y - X \beta\|_{\ell_2}^{2}  + \lambda \|\beta\|_{\ell_1} \right]$

Lassco stands for least absolute shrinkage and selection oberator, it performs both sgrinkage and soft-threshholding. The use of term $|\beta|_{\ell_1}=|\beta|$ is elegant in that, for $|\beta|^t$, when t>1 (like $|\beta|_{\ell_2}=|\beta|^2$) there is only change of slope in the primal form, t<1 it is not convex, only when t=1, there is change of slope and sharp change at cut point, while it is convex. Thus ee have $\hat{\beta}_{k,lassco}=max(0,{|\hat{\beta}_{k,ls}|}-\frac{\lambda}{|x_k^2|})$. 

In computation implementation, we use coordinate descent method. That is, start from 0 to iteratively update $\beta_k$ one at a time with $\Delta \beta$, only those represents large enough correlation will be admitted. The trick is to iterate over a list of $\lambda$ from large to small, so as $\lambda$ gets smaller more $\beta_k$ gets admitted. Thus the solution path of $\beta_k$ could show clearly how significant is each paramter. Lassco is suitable for sparse data.

## Classification

### Logistic Regression 

In Logistic Regression we deal with the circumstance where the response variable is not continuous, which is quite common in classification problems. Firstly we use a logit function to project it to a continous space, which can be interpretated as a score, and then the classification of an observation will depend on this score through the inverse of logit function --- sigmoid function:
\begin{align*}
{logit}(p_i) = \log(\frac{p_i}{1-p_i}) = X_i^T \beta \\  
\sigma(X_i^T \beta) =\frac{1}{1 + e^{-X_i^T \beta}}= \frac{e^{X_i^T\beta}}{1 + e^{X_i^T\beta}}
\end{align*}

With maximum likelihood estimation, we try to maximize likelihood function (often take log for simplicity in computation $l(\beta)$). We can derive two algorithm: gradient ascent and Newton-Ralphson, the latter is preferred due to its efficiency. 

#### IRLS (N-R method, gradient ascent) [[*View Code*]](#logistic)
Gradient ascent is to maximize $l(\beta)$) by find interatively updating $\beta$ with the product of lerning rate and the gradient (1st order derivative), learning rate governs the step size of each iteration. The Newton-Ralphson method improves this algorithm by using a more flexible learning rate $l''(\beta)^{-1}$). It is the curvature of $l$ at the point, and the step size should be smaller if it's large. With the weight $w_i=p_i(1-p_i)$, we can update $\beta$ as follows:   
\begin{align*}
\beta^{(t+1)} & = \left( \sum_{i=1}^{n} w_i X_iX_i^T \right)^{-1} ( \sum_{i=1}^{n} w_i X_i (X_i^T \beta^{(t)} + \dfrac{y_i-p_i}{w_i}))
\end{align*}

### Perceptron with +/- outcomes
In the context of training for classification, we use $y_i \in {+1,-1}$ instead of $y_i \in {+1,0}$, and the goal is to predict the classifier $\beta$. We start from $\beta_0=0$, and compute the score $X_i^T\beta$. If the score is positive, we classify y_i to be positive $(\hat{y_i}=sign(X_i^T\beta))$, vice versa. Then $y_iX_i^T\beta$ is the margin for ith observation, a negative margin indicates a mistake, the smaller it is the bigger the mistake. Therefore in the training process we learn from mistakes by minimize the total loss represented by the loss fucntion. This is called perceptron.

### Adaptive Boosting [[*View Code*]](#adaboost)
Adaboost means adaptive boosting. It is a committee machine, which consists of a number of base classifiers and the final classification is a perceptron based on the base classifiers. The training process miminze exponential loss: $y_i = {\rm sign}\left(\sum_{k=1}^{d} \beta_k h_k(X_i)\right)$,
where $\beta_k$ can be interpreted as the weight of vote of classifier $k$, and the goal is to derive a set of weights for optimal classfication result.

The traning employs coordinate descent, to sequentially add members to the committee, the weight is determined by how much error $h_\textrm{new}$ made on the weighted data. As a result, Adaboot could do on-line training with up-to-time data. 

### Support Vector Machine [[*View Code*]](#svm)
SVM can be interprated from the geometry perspective. Given the positive and negative examples in a high-dimension space spanned by all known features, the two part can beseparated by many separating hyperplanes. The goal is to choose the one with the maximum margin in order to guard against the random fluctuations in the unseen testing examples. Hence, we use a hinge loss $Hinge \; loss(\beta) = \sum_{i=1}^{n} max(0, 1-y_i X_i^T\beta)$. It has the advantage that, it penalizes not only the mistakes (negative margin), but also uncertainty (positive margin that is not large enough). Or say, it uses more information from the feature instead of deriving a binary weak classifier. In geometric view, it can separate the emaples with any hyperplanes while adaboost could only access to hyperplanes parallel to the axises.

The traning employs gradient descent.

### Neural Network [[*View Code*]](#nn)
A perceptron seeks to separate the positive examples and negative examples linearly, when they are not linearly separable, We may need to transform the original variables into some features so that they can be linearly separated. 

A way of understaning is to linearly map the original features to 1st layer hidden variables, and then to the response variable. In dealing with nonlinear approximation, neural net can also be viewed as a high-dimensional spline. If there are many layers in the neural net, the number of linear pieces is exponential in the number of layers. It can approximate highly non-linear mapping by patching up the large number of linear pieces. Therefore it can be vey powerful with large amount of training data. 

Neural Network is a generalized multi-layer perceptrons, logistic regression on top of logistic regressions. The gradient is then calculated by chain rule, so the chain back-propagates the error to assign the blame back to parameters on each level to updat eestimators.

# Tools

## R

### R and Rstudio
R is a widely used open source programming language and software environment for statistical computing and graphics. It is easy to use for starters, while capable of advanced programming supporting procedural programming with functions and object-oriented programming with generic functions. Besides, there is a vibrate on-line community with various open-source packages. These packages facilitate specialized statistical techiques including model-fitting, graphical device and reporting, and also provides good examples for coding. 

To get help in R, use "help" function. To improve R coding techniques, one can read source codes with "page" function and "getAnywhere" function. There are also various resources for self-learning:     
R news and tutorials contributed by (750) R bloggers: <https://www.r-bloggers.com/>.     
Developer community to find solutions to questions in programming: <https://stackoverflow.com/>.   

Rstudio is an interactive environment for R. It provides a user-friendly interface, and makes coding and editing in R and interaction with other languages a lot easier. It could directly compile Rcpp files, supports creation of R package and interactive web applications, easily connects to github and also provide a powerful tool for writing --- Rmarkdown.    
This page sumarizes keyboard shorcuts in Rstudio:<https://support.rstudio.com/hc/en-us/articles/200711853-Keyboard-Shortcuts>.     
```{r, echo=F, eval=F}
# Get help in R
? functionname

# View codes of existing functions
functionname
page(functionname)
## class function 
method(plot)
## function with a * or .
getAnywhere(plot.lm)
getS3method("print","lm")
getAnywhere(.login)
```

### R parallel
A major drawback of R is that it is not a fast language compared to other programming languages. It is because that R is both a language and an implementation of that language, it was designed to simplify coding for statistical analysis, but not to make things easier for computer.   

However, there are ways to make work more efficient by improving coding:   

1. Use and build on simple functions, specify the type and think more mathematically in array-oriented way.
```{r, echo=F, eval=F}
# instead of
colSums(object)
# use to save the configuring time
matrix(rep(1,n), nrow=1) %*% object
```

2. Use Grouping functions: tapply, by, aggregate, *apply family. It avoids long loop and makes codes shorter and more readable. This page is a helpful discussion on stackflow on these function: <https://stackoverflow.com/questions/3505701/grouping-functions-tapply-by-aggregate-and-the-apply-family>. Also, the package dpyr is a powerful tool in data cleaning with efficient backend: <https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html>

3. R parallel Package
```{r, echo=F, eval=F}
# package
library(doParallel)
# Setup parallel backend to use multiple processors
cl <- makeCluster(4)
registerDoParallel(cl)
# Start program
func <- foreach(i = 1:B, .combine = cbind) %dopar% {
       # combine results
       results
}
stopCluster(cl)
## it follows a similar logis to function lapply

## bootstrap example from homework
bootstrap_CI <- function(X, B = 10000){
  n <- length(X)
  ## Setup parallel backend to use multiple processors
  cl <- makeCluster(4)
  registerDoParallel(cl)
  ## Start program
  boostrap_results <- foreach(i = 1:B, .combine = cbind) %dopar% {
    boot_indices  <- sample(n,n,replace = T)
    boot_sample   <- X[boot_indices]
    boot_stat     <- median(boot_sample)
    result <- list(median = boot_stat)
    result
  }
  stopCluster(cl)
  my_results <- as.numeric(unlist(boostrap_results))
  CI_low     <- quantile(my_results,probs = 0.025)
  CI_high    <- quantile(my_results,probs = 0.975)
  CI         <- matrix(c(CI_low,CI_high),nrow = 2)
  return(CI)
}
```

Use package "doParallel" to setup parallel backend to use multiple processors and then start a program. The coding follows a similar logis to function "lapply". Some more references: <https://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/>,<https://cran.r-project.org/web/views/HighPerformanceComputing.html>.

### Rcpp
Another way to improve efficiency is to integrate C++ Code with package Rcpp. C++ is a fast language especially good at doing loop.        
Creat a cpp file in Rstudio, begin with the following and write a function in C++.
```{r, echo=F, eval=F}
#include <Rcpp.h>
using namespace Rcpp;

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
```
Source the file and then use the function in the environment. Some examples are included in the thrid part of this memo. For use of Armadillo, refer to page: <http://arma.sourceforge.net/docs.html#randu_randn_member>.

### R Package
By using R package, one can wrap all self-written utility functions scattered across numerous files. It allows a clear documentation of function, helps to use these function whenever it is needed and facilitate communications with partners. This page provides a detailed instruction with referneces to more advanced content: <http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html>.

Use package devtools could creat a package easily, package roxygen2 helps edit a standard ducumentation (remember to include #' @export!). Things can be a bit unstable in Rstudio when making a package, simply restarting the session may solve the problem. 
```{r, echo=F, eval=F}
# import packages
library(roxygen2)
library(devtools)
# creat the project
setwd("/Users/alice/Documents/R/StatsProgramming_202A")
create("LMjw")
# if the package requires rcpp
library(Rcpp)
library(RcppArmadillo)
RcppArmadillo.package.skeleton("JWlm",force=T)
# edit and save codes
## include #' @export!
# export functions and generate documents
setwd("/Users/alice/Documents/R/StatsProgramming_202A/LMjw")
document() 
compileAttributes(verbose=TRUE)
# install and try
setwd("/Users/alice/Documents/R/StatsProgramming_202A")
install("LMjw")
## or use Build-more-source package to create a .tar file
```

## Python
Python is a popular programming language, it is for general-purpose programming unlike R (for statistics use mainly). Python has a design philosophy that emphasizes code readability, and a syntax that allows fewer lines than that in languages such as C++ or Java. Python features a dynamic type system and automatic memory management, and supports multiple programming paradigms, including object-oriented, imperative, functional and procedural, and has a large and comprehensive standard library. This aunswer on Quora gives a neutral comparison between Python and R: <https://www.quora.com/Which-is-better-for-data-analysis-R-or-Python-Is-R-still-a-better-data-analysis-language-than-Python-Has-anyone-else-used-Python-with-Pandas-to-a-large-extent-in-data-analysis-projects>.

Data analysis in R and Python is actually quite similar, though R can be more convenient with built-in functions while Python relies on packages for mathematical use. Another difference is that index in Python starts from 0 instead of 1, which should be cautioned. However, Python is powerful for multiple uses and is sought after in tech industry. 

Anaconda is a open source distribution of the Python (and R) language for large-scale data processing. It incorporates Spider: a python development environment that works a lot like Rstudio; Jupiter: based on pytghon but runs code in many programming languages including R. For self-learning, stackflow is also useful, while there is an on-line course on edX for systematic study.

## SQL
SQL stands for Structured Query Language, is a language for managing data held in a relational database management system (RDBMS), or for stream processing in a relational data stream management system (RDSMS). It is a declarative language easy to learn, and widely used in today's industry. 

For basic operations, refer to tutorial: <https://sqlzoo.net/>.     
A foundation in DBS, relational algebra is useful, refer to on-line course: <https://lagunita.stanford.edu/courses/DB/2014/SelfPaced/about>. 

## Unix/Linux 
UNIX is an operating system. By operating system, we mean the suite of programs which make the computer work. It is a stable, multi-user, multi-tasking system for servers, desktops and laptops. The most popular varieties of UNIX are Sun Solaris, GNU/Linux, and MacOS X. In MacOS X, one can directly open terminal to operate on the system. Knowledge of UNIX is required for operations without GUI.

For basic commands, refers to tutorial: <http://www.ee.surrey.ac.uk/Teaching/Unix/unixintro.html>. It is funndamental to understand the concept of working directory and files, then to switch between directories, view and move files. Then with knowledge of input/output, one can use pipes to conveniently edit files. For specific command one can always refer to help command and google.

## Github
GitHub is a Web-based Git version control repository hosting service. It provides access control, task management, collaboration features to better manage personal or collaborative projects. With the large amount of open-source codes it hosts, it is a great help for programming study.

Each project on Github is saved as an on-line repository, that the user or users can pull it to local device, make and commit modification, and push modified files to the on-line repository. Guthub also allows comparison between different versions, and multiple users could edit on different sub-branch of a repository before merge it to the master branch.

### Use from terminal
GitHub can be used from terminal in unix.

To set up, one can generat a new SSH key and add it to the ssh-agent so that this device is linked to the github account, hence there is no need to reenter passphrase every time: 
<https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/>.

From the terminal, one can initialize a repository in a directory on the device, add files in this directory to the repository (maybe ignore some). Then with simple commands one can edit a file, compare it with previous version, add to repository, commit and push, also do branching and merges. This page offers a summaries for commands: <http://kbroman.org/github_tutorial/>.

### Use with Rstudio
Github can be used with Rstudio in straightforward steps following the instruction: <https://www.r-bloggers.com/rstudio-and-github/>. 

### Build a Website
Github can also host a websites for personal page or project page from a GitHub repository: <https://pages.github.com/>. One can basically edit and push everything from Rstudio with markdown files. The simplist way to start off is to fork from a friend's website<>, change the configuration, and add contents and features brick by brick. <https://github.com/Alice86/alice86.github.io/>

## Hoffman2
Hoffman2 is the UCLA shared cluster for huge data and/or heavy-duty computation. When log onto Hoffman2, one is assigned a "login node", and allowed access to 20GB storage in home directory. One can also write to a 2TB 7-day temporary storage /u/scratch, or use group space purchased by the project team.

To request an interactive session with specificed memory and time using the ‘qrsh’ command or submit a job on he schedule, a "compter node" is automatically assigned for execution. You will be informed via e-mail after the job is done or it is terminated.
```{r, echo=T, eval=F}
# Request a interactive session, requires 2 CPU cores, 4 GB of memory, 8 hours run time
qrsh -l h_rt=8:00:00,h_data=2G -pe shared 2
# Submit a job, requires 64 MB of memory, 8 hours run time, 
## named ‘gene_data_lasso’, refer to the job as ‘myjob.sh’.
## change the working directory to where you currently are in the file system -cwd
## emailed when the job begins-b, ends-e, or aborted-a
qsub -cwd -N gene_data_lasso -l h_data=64M,h_rt=08:00:00 -m bea myjob.sh
```

## Tensorflow
TensorFlow is an open source software library for numerical computation using data flow graphs in Python. It is for the purposes of conducting machine learning and deep neural networks research, but the system is general enough to be applicable in a wide variety of other domains as well. For installation, refers to <https://www.tensorflow.org/install/>. 
For a simple example of 2 layer neural network, view codes in the third part. [[*View Code*]](#tf)

## SAS
SAS is a software suite for statistical modeling, advanced analytics and data management, commonly used in business. It can be used with minimal programming knowledge, as it even has a point-and-click interface, and it generates well-formated results for reporting and sharing.
```{sas, echo=F, eval=F}
* Compute the correlation
proc corr data=data;
	var length mpg;
run;
* Make a scatterplot
proc sgplot data=data;
	scatter x = price y = mpg;
run;
* Make a box plot
proc boxplot data=data;
	plot mpg*foreign;
run;
* Perform simple linear regression without intercept term;
proc reg data=data;
	model mpg = price1 / noint;
run;
* Perform Generalized linear regression (polynomial); 
proc glm data=data;
	model mpg = length length*length;
run;
```

# Codes

## OLS Regression

### Gauss Jordan Elimination {#GJ}
```{r, echo=T, eval=F}
# Gauss Jordan elimination
## on invertible matrix A with m steps (if m=ncol(A) the result is inverse A)
## returns a matrix with the identity matrix on the left
myGaussJordan <- function(A, m){
       n<-dim(A)[1]
       B<-cbind(A,diag(rep(1,n)))
       for (k in 1:m) {
              # a <- B[k,k]
              # for (j in 1:(n*2)){
              #        B[k,j] <- B[k,j]/a
              # } 
              ## vector version, more efficient
              B[k,] <- B[k,]/B[k,k] 
              for (i in 1:n) {
                     if (i != k) {
                            b <- B[i,k]
                            # for (i in 1:(n*2)) {
                            #        B[i,j] <- B[i,j]-b*B[k,j]
                            # } 
                            B[i,] <- B[i,]-b*B[k,]
                     }
              }
       }
       return(B)
}
# LS Regression powered by Gauss Jordan
myLM_GJ <- function(X, Y){
       n <- nrow(X)
       p <- ncol(X)
       Z <- cbind(rep(1,n),X,Y)
       A <- t(Z)%*%Z
       S <- myGaussJordan(A,p+1)
       beta_hat <- S[1:(p+1),p+2]
       return(beta_hat)
}
```

Python
```{python, echo=T, eval=F}
import numpy as np
def myGaussJordanVec(A, m):
  n = A.shape[0]
  C = np.hstack((A, np.identity(n)))
  for k in range(m):
     # a = B[k,k]
     # for j in range(n*2):
     # B[k,j] = B[k,j]/a
     C[k,:] = C[k,:]/C[k,k]
     for i in range(n):
        if i !=k:
            C[i,:]=C[i,:]-C[k,:]*C[i,k]
  return B
# liner regression  
def myLinearRegression(X, Y):
  n = X.shape[0]
  p = X.shape[1]                         
  Z = np.hstack((np.repeat(1,n).reshape(n,1),X,Y))
  A= np.dot(Z.T,Z)
  S= myGaussJordanVec(A,p+1)
  beta_hat = S[0:p+1,p+1]
  return beta_hat
```

### Sweep Operator {#swp}
```{r, echo=T, eval=F}
# Sweep operation 
## on square matrix A with m steps
mySweep <- function(A, m){
       B <- A
       n <- dim(B)[1]
       for  (k in 1:m) {
              for (j in 1:n) {
              for (i in 1:n) {
                     if (i != k & j != k) {
                            B[i,j] <- B[i,j]-B[i,k]*B[k,j]/B[k,k]
                            }
                     }}
              for (i in 1:n) {
                     if (i != k) {
                            B[i,k] <- B[i,k]/B[k,k]
                     }
              }
              for (j in 1:n) {
                     if (j!=k) {
                            B[k,j] <- B[k,j]/B[k,k]
                     }      
              } 
              B[k,k] <- - 1/B[k,k]
       }
       return(B)
}
# LS Regression powered by Sweep
myLM_SWP <- function(X, Y){
       n <- nrow(X)
       p <- ncol(X)
       Z <- cbind(rep(1,n),X,Y)
       A <- t(Z) %*% Z
       S <- mySweep(A,p+1)
       beta_hat <- S[1:(p+1),p+2]
       return(beta_hat)
}
```

Rcpp
```{r engine='Rcpp', eval=F}
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat mySweepC(const mat A, int m) {
      // Declare our output
       mat B = A;
       int n = B.n_rows;
       for (int k=0; k<m; k++) {
              for (int j=0; j<n; j++) {
                     for (int i=0; i<n; i++) {
                            if (i != k & j != k) {
                                   B(i,j) = B(i,j)-B(i,k)*B(k,j)/B(k,k);
                            }
                     }}
              for (int i=0; i<n; i++) {
                     if (i != k) {
                            B(i,k) = B(i,k)/B(k,k);
                     }
              }
              for (int j=0; j<n; j++) {
                     if (j!=k) {
                            B(k,j) = B(k,j)/B(k,k);
                     }      
              } 
              B(k,k) = -1/B(k,k);
       }
       return(B);
}
// [[Rcpp::export()]]
mat myLinearRegressionC(const mat X, const mat Y){
       X: an 'n row' by 'p column' matrix of input variables.
       Y: an 'n row' by '1 column' matrix of responses
       int n = X.n_rows;
       int p = X.n_cols;
       mat O;
       O.ones(n,1);
       mat M = join_rows(O,X);
       mat Z = join_rows(M,Y);
       mat A = Z.t() * Z;
       mat S = mySweepC(A,p+1);
       mat beta_hat = S.submat(span(0,p),span(p+1,p+1));
       return(beta_hat);
}
```

Python
```{python, echo=T, eval=F}
import numpy as np

def mySweep(A, m):
    B = np.copy(A)   
    n = B.shape[0]
    for k in range(m):
        for i in range(n):
            for j in range(n):
                if i<>k and j<>k:
                    B[i,j] = B[i,j]-B[i,k]*B[k,j]/B[k,k];                
        for i in range(n):
            if i<>k:
                B[i,k]=B[i,k]/B[k,k];
        for j in range(n):
            if j<>k:
                B[k,j]= B[k,j]/B[k,k];
        B[k,k] = - 1/B[k,k];
    return(B)
```

### QR Decomposition {#qr}
```{r, echo=T, eval=F}
# QR Decomposition
## on symmetric matrix A
myQR <- function(A){
       n <- dim(A)[1]
       m <- dim(A)[2]
       R <- A
       Q <- diag(n)
       for (i in 1:(m-1)) {
              X <- matrix(rep(0,n),nrow = n)
              X[i:n] <- R[i:n,i]
              V <- X
              V[i] <- X[i]+norm(X,"F")*sign(X[i])
              U <- V/norm(V,"F")
              R <- R-2*(U%*%t(U)%*%R)
              Q <- Q-2*U%*%t(U)%*%Q
       }
       return(list("Q" = t(Q), "R" = R))
}
# linear regression
myLM <- function(X, Y){
       n <- nrow(X)
       p <- ncol(X)
       A <- cbind(rep(1,n),X,Y)
       
       R <- myQR(A)$R
       R1 = R[1:(p+1),1:(p+1)]
       Y1 = R[1:(p+1),p+2]
       beta_ls <- solve(R1)%*%Y1
       return(beta_ls,sigm.sq)
}
```

Rcpp
```{r engine='Rcpp', eval=F}
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
double signC(double d){
       return d<0?-1:d>0? 1:0;
}

// [[Rcpp::export()]]
List myQRC(const mat A){ 
       int n = A.n_rows;
       int m = A.n_cols;
       mat R = A;
       mat Q = eye<mat>(n,n);
       for (int i=0; i<m-1; i++) {
              mat X;
              X.zeros(n,1);
              X(span(i,n-1),0) = R(span(i,n-1),i);
              mat V = X;
              double x = X(i,0);
              double sign = signC(x);
              double nrm = norm(X);
              V(i,0) = x+sign*nrm;
              double nrmv = norm(V);
              mat U = V/nrmv;
              R -= 2 * U * U.t() * R;
              Q -= 2 * U * U.t() * Q;
       }
       output["Q"] = Q.t();
       output["R"] = R;
       return(output);
}
# least square
// [[Rcpp::export()]]
mat myLinearRegressionC(const mat X, const mat Y){
       int n = X.n_rows;
       int p = X.n_cols;
       mat O;
       O.ones(n,1);
       mat M = join_rows(O,X);
       mat A = join_rows(M,Y);
       
       mat R = myQRC(A)["R"];
       mat R1 = R(span(0,p),span(0,p));
       mat Y1 = R(span(0,p),(p+1));
       mat beta_ls = (R1.i() * Y1).t();
       return(beta_ls.t());
} 
```

Python
```{r engine='Rcpp', eval=F}
import numpy as np 
from scipy import linalg
# QR Decomposition
def qr(A):
    n,m=A.shape
    R=A.copy()
    Q=np.eye(n)
    for k in range(m-1):
        X = np.zeros((n,1))
        X[k:]=R[k:,k]
        V = X
        V[k]=X[k]+np.sign(X[k])*np.linalg.norm(X)
        S=np.linalg.norm(V)
        U=V/S
        R -= 2*np.dot(U,np.dot(U.T,R))
        Q -= 2*np.dot(U,np.dot(U.T,Q))
    Q=Q.T
    return Q,R
```

### PCA (eigen decomposition){#pca}
```{r, echo=T, eval=F}
# Eigen decomposition
## on symmetric matrix A, interation times for power method numIter default to 1000
## return eigen values D and vectors V
myEigen_QR <- function(A, numIter = 1000){
       r <- dim(A)[1]
       v <- matrix(rnorm(r^2),nrow=r)
       for (i in 1:numIter) {
              QR <- myQR(v)
              Q = QR$Q
              R = QR$R
              v = A %*% Q
       }
       return(list("D" = diag(R), "V" = Q))
}
```

Rcpp
```{r engine='Rcpp', eval=F}
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
List myEigen_QRC(const mat A, const int numIter = 1000){
       mat B = A;
       int p = A.n_cols;
       mat V;
       V.randn(p,p);
       for (int i=0; i<=numIter; i++) {
              mat Q = myQRC(V)["Q"] ;
              V = B * Q;
       }
       mat Q = myQRC(V)["Q"] ;
       mat R = myQRC(V)["R"] ;
       vec D = R.diag();

       output["D"] = D;
       output["V"] = Q;
       return(output);
}
```

Python
```{python, echo=T, eval=F}
import numpy as np 
# eigen decomposition
def eigen_qr(A):
    T = 1000
    A_copy = A.copy
    r,c = A_copy.shape
    v = np.random.random_sample((r,r))
    for i in range(T):
        Q,R = qr(v)
        v =np.dot(A_copy,Q)
 
    return R.diagonal(),Q
```

## Regularized Learning

Test Data
```{r, echo=T, eval=F}
import numpy as np
import sklearn.datasets as ds
from sklearn.model_selection import train_test_split
def prepare_data(valid_digits=np.array((6, 5))):
    if len(valid_digits) != 2:
        raise Exception("Error: you must specify exactly 2 digits for classification!")
    data = ds.load_digits()
    labels = data['target']
    features = data['data']
    X = features[(labels == valid_digits[0]) | (labels == valid_digits[1]), :]
    Y = labels[(labels == valid_digits[0]) | (labels == valid_digits[1]),]
    X = np.asarray(map(lambda k: X[k, :] / X[k, :].max(), range(0, len(X))))
    Y[Y == valid_digits[0]] = 0
    Y[Y == valid_digits[1]] = 1
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.25, random_state=10)
    Y_train = Y_train.reshape((len(Y_train), 1))
    Y_test = Y_test.reshape((len(Y_test), 1))
    return X_train, Y_train, X_test, Y_test
```

### Ridge Regression {#ridge}
```{r, echo=T, eval=F}
myRidge <- function(X, Y, lambda){
       n <- nrow(X)
       p <- ncol(X)
       Z <- cbind(rep(1, n), X, Y) 
       A <- t(Z) %*% Z
       D <- diag(rep(lambda, p+2)) 
       D[1, 1] <- 0
       D[p+2, p+2] <- 0
       A <- A+D
       S <- mySweep(A, p+1)
       beta_ridge <- S[1:(p+1), p+2]
       return(beta_ridge)
}
```

### Spline Regression {#spline}
```{r, echo=T, eval=F}
myRidge <- function(X, Y, lambda){
       n <- nrow(X)
       p <- ncol(X)
       Z <- cbind(rep(1, n), X, Y) 
       A <- t(Z) %*% Z
       D <- diag(rep(lambda, p+2)) 
       D[1, 1] <- 0
       D[p+2, p+2] <- 0
       A <- A+D
       S <- mySweep(A, p+1)
       beta_ridge <- S[1:(p+1), p+2]
       return(beta_ridge)
}
```

### Lasso Regression {#lassco}
```{r, echo=T, eval=F}
myLasso <- function(X, Y, lambda_all){
       n <- nrow(X)
       X <- cbind(rep(1,n),X)
       p <- ncol(X)
       L <- length(lambda_all)
       SS <- rep(0,p)
       for(j in 1:p) {
              SS[j] <- sum(X[ ,j]^2)
       }
       S <- SS
       S[1] <- Inf
       beta_all <- matrix(rep(0,p*L), nrow=p)
       beta <-  matrix(rep(0,p), nrow=p)
       R <- Y
       for(l in 1:L) {
              lambda <- lambda_all[l]
              for(t in 1:10) {
                     for(k in 1:p) {
                            db <- sum(R*X[,k])/SS[k]
                            b <- beta[k] + db
                            b <- sign(b)*max(0, abs(b)-lambda/S[k])
                            db <- b-beta[k]
                            R <- R - X[ ,k]*db
                            beta[k] <- b
                     }
              }
              beta_all[,l] <- beta
       }
       return(beta_all)
}
```

## Classification

### Logistic Regression {#logistic}
```{r, echo=T, eval=F}
## Expit/sigmoid function
expit <- function(x){
       1 / (1 + exp(-x))
}
# Logistic regression 
## IRLS powered by myLinearRegressionC
myLogistic <- function(X, Y){
       n <- nrow(X)
       m <- ncol(X)
       
       beta <- matrix(rep(0,m),nrow = m)
       epsilon <- 10^(-6)
       err <- 1
       while (err >= epsilon) {
              eta <- X %*% beta
              p <- expit(eta)
              w <- p*(1-p)
              Y_hat <- eta + (Y-p)/w
              
              X_ls <- matrix(rep(sqrt(w),each=m),byrow=T,ncol = m) * X
              Y_ls <- sqrt(w) * Y_hat
              beta_n <- myLM(X_ls, Y_ls)
              err <- sum(abs(beta_n-beta))
              beta <- beta_n
       }
       beta  
}
```

### Adaptive Boosting {#adaboost}
```{r, echo=T, eval=F}
# Adaptive Boosting
## weak classifier: X_train > 0.8
myAdaboost <- function(X_train, Y_train, X_test, Y_test,
                       num_iterations = 200) {
  n <- dim(X_train)[1] 
  p <- dim(X_train)[2]
  threshold <- 0.8
  X_train1 <- 2 * (X_train > threshold) - 1
  Y_train  <- 2 * Y_train - 1
  X_test1 <- 2 * (X_test > threshold) - 1
  Y_test  <- 2 * Y_test - 1
  beta <- matrix(rep(0,p), nrow = p)
  w <- matrix(rep(1/n, n), nrow = n)
  weak_results <- Y_train * X_train1 > 0
  acc_train <- rep(0, num_iterations)
  acc_test  <- rep(0, num_iterations)
  for(it in 1:num_iterations) {
    w <- w / sum(w)
    weighted_weak_results <- w[,1] * weak_results
    weighted_accuracy <- colSums(weighted_weak_results)
    e <- 1 - weighted_accuracy
    j <- which.min(e)
    dbeta <-log((1-e[j])/e[j])/2
    beta[j] <- beta[j] + dbeta
    w <-  w*exp(-y*x[, j]*dbeta)
    acc_train[it] <- mean(sign(X_train1 %*% beta) == Y_train)
    acc_test[it] <- mean(sign(X_test1 %*% beta) == Y_test)   
  }
  output <- list(beta = beta, acc_train = acc_train, acc_test = acc_test)
  output
}
```

Python
```{python, echo=T, eval=F}
def my_Adaboost(X_train, Y_train, X_test, Y_test, num_iterations=200):
    n = X_train.shape[0]
    p = X_train.shape[1]
    threshold = 0.8
    X_train1 = 2 * (X_train > threshold) - 1
    Y_train = 2 * Y_train - 1
    X_test1 = 2 * (X_test > threshold) - 1
    Y_test = 2 * Y_test - 1
    beta = np.repeat(0., p).reshape((p, 1))
    w = np.repeat(1. / n, n).reshape((n, 1))
    weak_results = np.multiply(Y_train, X_train1) > 0

    acc_train = np.repeat(0., num_iterations, axis=0)
    acc_test = np.repeat(0., num_iterations, axis=0)
    for it in range(0,num_iterations):
        w = w/sum(w)
        a = np.dot(np.repeat(1,n).reshape(1,n),w * weak_results)
        e = 1-a
        k = np.argmin(e)
        db = np.log((1-e[0,k])/e[0,k])/2
        beta[k,0] = beta[k,0]+db
        w = (w[:,0] * np.exp(-Y_train[:,0]*X_train1[:,k]*db)).reshape(n,1)
        acc_train[it] = np.mean(np.sign(np.dot(X_train1, beta)) == Y_train)
        acc_test[it] = np.mean(np.sign(np.dot(X_test1, beta)) == Y_test)
    return beta, acc_train, acc_test
```

### Support Vetor Machine {#svm}
```{r, echo=T, eval=F}
# SVM
my_SVM <- function(X_train, Y_train, X_test, Y_test, lambda = 0.01,
                   num_iterations = 1000, learning_rate = 0.1){
  n <- dim(X_train)[1]
  p <- dim(X_train)[2] + 1
  X_train1  <- cbind(rep(1, n), X_train)
  Y_train <- 2 * Y_train - 1
  beta <- matrix(rep(0, p), nrow = p)
  ntest <- nrow(X_test)
  X_test1 <- cbind(rep(1, ntest), X_test)
  Y_test <- 2 * Y_test - 1
  acc_train <- rep(0, num_iterations) 
  acc_test  <- rep(0, num_iterations) 
  for(it in 1:num_iterations) {
    s <- X_train1 %*% beta
    db <- s * Y_train < 1
    dbeta <- matrix(rep(1, n), nrow = 1) %*%((matrix(db*Y, n, p)*X1))/n; 
    beta <- beta + learning_rate * t(dbeta)
    beta[2:p] <- beta[2:p] - lambda * beta[2:p]
    acc_train[it] <- mean(sign(X_train1 %*% beta) == Y_train)
    acc_test[it] <- mean(sign(X_test1 %*% beta) == Y_test)
    }
  model <- list(beta = beta, acc_train = acc_train, acc_test = acc_test)
  model
}
```

Python
```{python, echo=T, eval=F}
import numpy as np
def my_SVM(X_train, Y_train, X_test, Y_test, lamb=0.01, num_iterations=200, learning_rate=0.1):
    n = X_train.shape[0]
    p = X_train.shape[1] + 1
    X_train1 = np.concatenate((np.repeat(1, n, axis=0).reshape((n, 1)), X_train), axis=1)
    Y_train = 2 * Y_train - 1
    beta = np.repeat(0., p, axis=0).reshape((p, 1))
    ntest = X_test.shape[0]
    X_test1 = np.concatenate((np.repeat(1, ntest, axis=0).reshape((ntest, 1)), X_test), axis=1)
    Y_test = 2 * Y_test - 1
    acc_train = np.repeat(0., num_iterations, axis=0)
    acc_test = np.repeat(0., num_iterations, axis=0)
    for it in range(0,num_iterations):
        score = np.dot(X_train1, beta)
        delta = score * Y_train < 1
        dbeta = np.sum(delta * np.multiply(Y_train, X_train1),axis=0)/n
        beta = beta + learning_rate * dbeta.reshape(p,1)
        beta[1:(p-1),0] = beta[1:(p-1),0] - lamb * beta[1:(p-1),0]
        acc_train[it] = np.mean(np.sign(np.dot(X_train1, beta)) == Y_train)
        acc_test[it] = np.mean(np.sign(np.dot(X_test1, beta)) == Y_test)
    return beta, acc_train, acc_test
```

### Neural Network{#nn}
```{r, echo=T, eval=F}
# Neoral Networl
## 2 layers with sigmoid ativation function
my_NN <- function(X_train, Y_train, X_test, Y_test, num_hidden = 20, 
                     num_iterations = 1000, learning_rate = 1e-1) {
       n <- dim(X_train)[1]
       p <- dim(X_train)[2] + 1
       ntest <- dim(X_test)[1]
       X_train1 <- cbind(rep(1, n), X_train)
       X_test1 <- cbind(rep(1, ntest), X_test)
       alpha <- matrix(rnorm(p * num_hidden), nrow = p)
       beta  <- matrix(rnorm((num_hidden + 1)), nrow = num_hidden + 1)
  
       acc_train     <- rep(0, num_iterations)
       acc_test      <- rep(0, num_iterations)
  
       for(it in 1:num_iterations) {
              Z  <- 1 / (1 + exp(-X_train1 %*% alpha))
              Z1 <- cbind(rep(1, n), Z)
              pr <- 1 / (1 + exp(-Z1 %*% beta))
              dbeta = matrix(rep(1, n), nrow = 1) %*%((matrix(Y_train-pr, n, m+1)*Z1))/n; 
              beta  <- beta + learning_rate * t(dbeta)
              for(k in 1:num_hidden) {
                     da <- (Y_train - pr)*beta[k+1]*Z[, k]*(1-Z[, k])
                     dalpha <- matrix(rep(1, n), nrow = 1)%*%((matrix(da, n, p)*X_train1))/n
                     alpha[, k] <- alpha[, k] + learning_rate * t(dalpha)
              }    
       acc_train[it] <- accuracy(pr, Y_train)
       Ztest         <- 1/(1 + exp(-X_test1 %*% alpha))
       Ztest1        <- cbind(rep(1, ntest), Ztest)
       prtest        <- 1/(1 + exp(-Ztest1 %*% beta))
       acc_test[it]  <- accuracy(prtest, Y_test)
       cat("On iteration ", it, " the training accuracy is ", acc_train[it], 
        " and the testing accuracy is ", acc_test[it], sep = "")
       }
       model <- list(alpha = alpha, beta = beta, acc_train = acc_train, acc_test = acc_test)
       model
}
```

Python
```{r, echo=T, eval=F}
# Neoral Networl
## 2 layers with Relu for 1st layer (alpha) and sigmoid for 2nd layer (beta)
def my_NN(X_train,Y_train,X_test,Y_test,num_hidden=20,num_iterations=1000,learning_rate=1e-1):
    n=X_train.shape[0]
    p=X_train.shape[1]+1
    ntest=X_test.shape[0]
    X_train1= np.concatenate((np.repeat(1,n,axis=0).reshape((n,1)),X_train),axis=1)
    X_test1 = np.concatenate((np.repeat(1, ntest, axis=0).reshape((ntest, 1)), X_test), axis=1)
    alpha=np.random.standard_normal((p,num_hidden))
    beta=np.random.standard_normal((num_hidden+1,1))
    acc_train=np.repeat(0.,num_iterations)
    acc_test=np.repeat(0.,num_iterations)

    for it in range(0,num_iterations):
        h1 = np.maximum(np.dot(X_train1, alpha),0)
        h = np.concatenate((np.repeat(1,n,axis=0).reshape((n,1)),h1),axis=1)
        pr = 1/(1+np.exp(-np.dot(h,beta)))
        dbeta = np.mean((Y_train-pr).reshape(n,1)*h,axis=0)
        beta = beta + (learning_rate*dbeta).reshape(num_hidden+1,1)
        for k in range(0,num_hidden):
            da = (Y_train - pr)[:,0]*beta[k+1,0]*np.sign(h1[:,k])
            dalpha = np.mean(da.reshape(n,1)*X_train1,axis=0)
            alpha[:,k] = alpha[:,k] + learning_rate*dalpha
        acc_train[it] = accuracy(pr, Y_train)
        ht = 1/(1+np.exp(-np.dot(X_test1, alpha)))
        htest = np.concatenate((np.repeat(1,ntest,axis=0).reshape((ntest,1)),ht),axis=1)
        prtest = 1/(1+np.exp(-np.dot(htest,beta)))
        acc_test[it] = accuracy(prtest, Y_test)    
        
    return alpha,beta,acc_train,acc_test
```


Tensorflow {#tf}
```{r, echo=T, eval=F}
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import argparse
import sys
import os
os.environ["CUDA_VISIBLE_DEVICES"]=""
from tensorflow.examples.tutorials.mnist import input_data
import tensorflow as tf

FLAGS = None
def main(_):
  # Import data
  mnist = input_data.read_data_sets(FLAGS.data_dir, one_hot=True)
  num_hidden1 = 100
  num_hidden2 = 10
  x = tf.placeholder(tf.float32, [None, 784])
  W1 = tf.Variable(tf.random_uniform([784, num_hidden1],-0.01,0.01))
  b1 = tf.Variable(tf.random_normal_initializer()([num_hidden1]))
  W2 = tf.Variable(tf.random_uniform([num_hidden1, num_hidden2],-0.01,0.01))
  b2 = tf.Variable(tf.random_normal_initializer()([num_hidden2]))
  z = tf.nn.relu(tf.matmul(x, W1) + b1)
  y = tf.nn.softmax(tf.matmul(z, W2) + b2)
  
  # Define loss and optimizer
  y_ = tf.placeholder(tf.float32, [None, 10])
  cross_entropy = tf.reduce_mean(-tf.reduce_sum(y_ * tf.log(y), reduction_indices=[1]))
  train_step = tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)
  sess = tf.InteractiveSession()
  tf.global_variables_initializer().run()
  for _ in range(1000):
      batch_xs, batch_ys = mnist.train.next_batch(100)
      sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})
      
  # test trained model  
  correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
  accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
  print(sess.run(accuracy, feed_dict={x: mnist.test.images, y_: mnist.test.labels}))
  sess.close()
  
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--data_dir', type=str, default='/tmp/tensorflow/mnist/input_data',
                      help='Directory for storing input data')
  FLAGS, unparsed = parser.parse_known_args()
  tf.app.run(main=main, argv=[sys.argv[0]] + unparsed)
```  
