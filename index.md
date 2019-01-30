---
layout: page
title: D-Optimality with Level Balance Constraints
tagline: Optimizing multiply-constrained matrices for choice research
description:
---
*Andrew Armstrong, Alex Graber*

# Background: What is Choice Research?

Pasta Example.  What is choice research?  

# Business Problem

Adding constraints

<details><summary>D-Optimality</summary>
  <div markdown = "1">
    
# D-Optimality
### Theory

### Fedorov Algorithm
Lay out base optimization problem

  </div>
</details>

<details><summary>Adding Constraints, Pt 1</summary>
  <div markdown = "1">
    
# Adding Constraints, Pt 1
### Toy Problem 


### Model Definition

# Constrained D-Optimality
### Theory

### Modified Fedorov Algorithm

### Parallelized

#### MFA Results

### Genetic Algorithm

#### GA Results

  </div>
</details>

<details><summary>Adding Constraints, Pt 2</summary>
  <div markdown = "1">
    
# Adding Constraints, Pt 2
### Full Problem

### Model Definition

### Results & Performance

  </div>
</details>

<details><summary>Bibliography</summary>
  <div markdown = "1">
    
# Bibliography

1. Hauser, J., & Rao, V. (2002, September). Conjoint Analysis, Related Modeling, and Applications. In IN MARKET RESEARCH AND MODELING: PROGRESS AND PROSPECTS: A TRIBUTE. Kluwer Academic Publishers.
2.  Kuhfeld, W., Huber, J., & Zwerina, K. (1996, September). A General Method for Constructing Efficient Choice Designs. Retrieved October 2018, from https://faculty.fuqua.duke.edu/~jch8/bio/Papers/Zwerina%20Kuhfeld%20Huber.pdf
3. Labadi, L. A. (2015, February). Some Refinements on Fedorov’s Algorithms for Constructing D-optimal Designs. Brazilian Journal of Probability and Statistics, 29, 53-70.
4. Triefenback, F. (2008). Design of Experiments: The D-Optimal Approach and Its Implementation As a Computer Algorithm. Umeå University, Department of Computing Science.
5. Warren F. Kuhfeld. (2001, January). Multinomial Logit, Discrete Choice Modeling. Retrieved October 2018, from https://www.stat.auckland.ac.nz/~balemi/Choice.pdf

  </div>
</details>


$$
\left[ \begin{array} { c c c } { A _ { 1 } } & { G _ { 1 } } & { B _ { 1 } } \\ { \vdots } & { \vdots } & { \vdots } \\ { A _ { N } } & { G _ { N } } & { B _ { N } } \end{array} \right]
$$

$$ A _ { i } , \quad A \in \{ 0,1,2 \} , i = 1 , \ldots , n $$

inline $$ 
\text { maximize } f ( X ) = 100 \frac { \operatorname { Det } \left( X ^ { T } X \right) ^ { 1 / p } } { N } - \lambda \sum \left| \delta _ { \text { distribution } } \right| - \lambda ^ { 2 } \sum \left| \delta _ { \text {interaction} } \right|
$$ inline


block
\\[
\text { maximize } f ( X ) = 100 \frac { \operatorname { Det } \left( X ^ { T } X \right) ^ { 1 / p } } { N } - \lambda \sum \left| \delta _ { \text { distribution } } \right| - \lambda ^ { 2 } \sum \left| \delta _ { \text {interaction} } \right|
\\]
block


{% raw %}
$$a^2 + b^2 = c^2$$ --> note that all equations between these tags will not need escaping! 
{% endraw %}

