---
layout: page
title: D-Optimality with Level Balance Constraints
tagline: Optimizing multiply-constrained matrices for choice research
description:
---
*Andrew Armstrong, Alex Graber*

## Background: What is Choice Research?

Let's say you're having trouble deciding what to have for dinner.  You know you want a pasta dish, but you don't know which dish - there are so many options!  You have to choose the noodle shape, the sauce, and what kind of topping you want.  


![Pasta](/assets/Picture1.png)  
*Pasta Options (source: google image search)*  


Choice research could help you make your decision.  Using the conjoint method of choice research, you could ennumerate the various options for each attribute and come up with a map defining a pasta dish based on the various attributes:


![Pasta Design](/assets/Picture2.png)  
*Pasta Choice Design*  


Notice how the pasta dish in the picture is dish #2 in the chart on the right; it is composed of Spaghetti (pasta level 0), Pesto (sauce level 0), and topped with Parmasen Cheese (topping option 1).  


If you had to try every possible pasta dish in the world, it would be impossible - there are so many types!  But by defining pasta by its attributes, you can select a small subset of the entire variety and use a multinomial logit model to understand which attributes are your favorite.  Then, you could make yourself your personal ultimate pasta dish by combining all of your favorite attributes.  




<details><summary>D-Optimality</summary>
  <div markdown = "1">

## D-Optimality
The pasta story is a simplistic example of choice research, but it should give you the intuition for why choice research is important, and how you can use a smaller portion of all possible options to associate value or importance with attribute levels.  This raises a key question: *How can you identify the best subset to use that maximizes the information gained from the research?*

It is clear that when the number of attributes and levels grow beyond a small set, presenting the full design (full factorial) becomes a challenge due to both the number of combinations required and the amount of burden placed on the respondent.  Fractional factorial designs, then, seek to allow the research to eke as much data out of the analysis as possible but use a much more limited subset of stimuli – but how do we know what the best (i.e., most efficient) fractional factorial design is?


Much research has been done on the topic of identifying efficient experimental designs (Hauser & Rao, 2002).  The current standard used to identify ‘efficient design’ is D-error – the geometric mean of the eigenvalues of the covariance matrix (D-efficiency is the inverse of D-error) (Kuhfeld, Huber, & Zwerina, 1996). Thus, the goal of an efficient design is to minimize D-error (therefore maximizing D-efficiency).  


D-efficient designs satisfy four principles (Kuhfeld, Huber, & Zwerina, 1996):
* **Orthogonality** is satisfied when the levels of each attribute vary independently of one another.  
* **Level balance** is satisfied when the levels of each attribute appear with equal frequency. 
* **Minimal overlap** is satisfied when the alternatives within each choice set have nonoverlapping attribute levels.  
* **Utility balance** is satisfied when the utilities of alternatives within choice sets are the same.


The standard method to identify an efficient design is to use one of any variant of the **Fedorov Algorithm** which, given a starting design, recursively makes exchange(s) that reduce D-error until some convergence criteria is met.  This method is susceptible to local minima; it may be necessary to run multiple iterations of the Fedorov Algorithm with different random starting designs to find the most efficient design (Kuhfeld, Huber, & Zwerina, 1996).

  </div>
</details>




## Business Problem

A pharmaceutical market research firm uses simulated patient treatment as a method to understand physician demand in specific treatment areas.  In this method, a limited universe of patients is designed in order to represent as much of the actual disease area patient universe as possible.  Patients are defined by multiple attributes (age, gender, BMI, etc.), and each attribute may have multiple levels (male/female, etc.).  The distribution of attribute levels is controlled such that the demographics of the simulated patient universe approximate the real patient universe.  These simulated patients are then treated, where a given treatment (yes/no) can be related back to the patient design.


Patient simulation research in this manner is a specialized choice methodology somewhat analogous to conjoint.  In both conjoint and patient simulation, respondents are forced to make a decision based on a stimulus that is composed of multiple attributes and levels (Kuhfeld, Huber, & Zwerina, 1996).  There are two key differences between conjoint (the pasta example) and patient simulation:
1. The design for patient simulation inherently contains *d-error* as a result of violating the principle of **level balance**:  Since the goal is for the simulated patient universe to map to the actual patient universe, the researcher may need to control for the distribution of levels within each attribute.  
2. Certain attributes and levels may have required interactions (i.e., a patient must be female to be pregnant), potentially violating **orthogonality** and **minimal overlap**.  




<details><summary>Adding Constraints, Pt 1</summary>
  <div markdown = "1">
    
## Adding Constraints, Pt 1
### Toy Problem 

As a toy problem, let us consider a patient universe in which patients are defined by:
![Picture 3](/assets/Picture3.png)


Expanding out all possibilities into the entire candidate set, this would be 3\*3\*2\*3 = 54 unique patient profiles.  Given respondent time is expensive, and high respondent burden decreases quality of results, we seek to reduce time-in-survey by creating a fractional-factorial design of 8 unique patient profiles.  As we want to extract as much data from the exercise as possible, the 8-profile fractional-factorial design must be as efficient as possible.


Practically speaking, the number of attributes is limited to no more than 25, each with at most 5 levels due to the complexity of the simulation, limited respondent pool, and limited number of experiments possible per respondent.  Thus, at most, the candidate set contains $$ 5 ^ { 25 } $$ (approx. $$ 3 \times 10 ^ { 17 } $$) possibilities – and will generally be significantly smaller as not all 25 attributes are used and most contain fewer than 5 levels.  However, the worst-case scenario requires approximately $$ 2 \times 10 ^ { 10 } $$ gigabytes to merely store the candidate set.  The combinatorics problem explains why stochastic search algorithms such as simulated annealing or genetic algorithms are frequently used instead of an exhaustive search against a complete candidate set.  


### Model Definition

Our goal is to maximize the weighted d-optimality of the design matrix, penalized for missing distributions and impossible variable interactions, and subject to the distributions of each attribute’s levels and interactions, where each attribute’s level is represented by a binary variable.

Objective Function (Wanida Limmun, 2012): 

$$ 
\text { maximize } f ( X ) = 100 \frac { \operatorname { Det } \left( X ^ { T } X \right) ^ { 1 / p } } { N } - \lambda \sum \left| \delta _ { \text { distribution } } \right| - \lambda ^ { 2 } \sum \left| \delta _ { \text {interaction} } \right|
$$ 

where $$ N $$ is the number of observations, $$ \delta $$ are vectors of relaxation variables, and $$ X $$ is the design matrix:

$$  
\left[ \begin{array} { c c c } { A _ { 1 } } & { G _ { 1 } } & { B _ { 1 } } \\ { \vdots } & { \vdots } & { \vdots } \\ { A _ { N } } & { G _ { N } } & { B _ { N } } \end{array} \right]
$$  

With decision variables $$ A, B, G $$ representing attributes:  

$$\qquad$$ Age:  $$ A _ { i } , \quad A \in \{ 0,1,2 \} , i = 1 , \ldots , n $$ *(i.e., the age group classification for each patient i)*
    
$$\qquad$$ Gender: $$ G _ { i } , \quad G \in \{ 0,1 \} , i = 1 , \ldots , n $$ *(i.e., the gender classification for each patient i)*
    
$$\qquad$$ BMI: $$ B _ { i } , \quad B \in \{ 0,1,2 \} , i = 1 , \ldots , n $$ *(i.e., the BMI classification for each patient i)*


For easier constraint formulation, we can use the Dantzig-Wolfe reformulation to rewrite our integer variables where the capital letter represents the binary variable series replacing an integer variable, and the lowercase letter represents the integer set of levels permissible for the given attribute:

$$  
\begin{array} { l l } { A _ { i } = \sum _ { 0 } ^ { z } z Z _ { z } , } & { \text { and } \sum _ { 0 _ { y } } ^ { z } Z _ { z } = 1 , \quad Z _ { z } \in \{ 0,1 \} , z \in \{ 0,1,2 \} } \\ { G _ { i } = \sum _ { 0 } ^ { y } y Y _ { y } , } & { \text { and } \sum _ { w ^ { 0 } } ^ { z } Y _ { y } = 1 , \quad Y _ { y } \in \{ 0,1 \} , y \in \{ 0,1 \} } \\ { B _ { i } = \sum _ { 0 } ^ { w } w W _ { w } , } & { \text { and } \sum _ { 0 } ^ { w ^ { 0 } } W _ { w } = 1 , \quad W _ { w } \in \{ 0,1 \} , w \in \{ 0,1,2 \} } \end{array}
$$  

Subject to:

$$\qquad$$ Age group proportions:

$$
\begin{array} { l } { \frac { \sum Z _ { 0 } } { N } = .25 + \delta _ { Z 0 } } \\ { \frac { \sum Z _ { 1 } } { N } = .5 + \delta _ { Z 1 } } \\ { \frac { \sum Z _ { 2 } } { N } = .25 + \delta _ { Z 1 } } \end{array}
$$

$$\qquad$$ Gender proportions:

$$
\begin{array} { l } { \frac { \sum Y _ { 0 } } { N } = .5 + \delta _ { Y 0 } } \\ { \frac { \sum Y _ { 1 } } { N } = .5 + \delta _ { Y 0 } } \end{array}
$$

$$\qquad$$ BMI proportions:

$$
\begin{array} { l } { \frac { \sum W _ { 0 } } { N } = .25 + \delta _ { W 0 } } \\ { \frac { \sum W _ { 1 } } { N } = .25 + \delta _ { W 1 } } \\ { \frac { \sum W _ { 2 } } { N } = .5 + \delta _ { W 2 } } \end{array}
$$

$$\qquad$$ Binary constraints: $$ W , Y , Z \in \{ 0,1 \} $$

$$\qquad$$ Interaction slacks: While not specified in the toy problem, it is entirely possible that we may have interactions specified in the design (i.e., men cannot be pregnant).  In these cases, the slacks are the count of the impossible interactions.  We will penalize these interaction slacks twice because they are more costly to the design than a missed distribution. 

### Constrained D-Optimality

For our discrete-choice design, the information matrix of an *n*-point design is $$\mathrm { M } = X _ { n } ^ { T } X _ { n } $$, where $$ X $$ is an $$ n \times p $$ design matrix. 

### Modified Fedorov Algorithm

### Parallelized

#### MFA Results

### Genetic Algorithm

#### GA Results

  </div>
</details>

<details><summary>Adding Constraints, Pt 2</summary>
  <div markdown = "1">
    
## Adding Constraints, Pt 2
### Full Problem

### Model Definition
  </div>
</details>



    
## Conclusion
### Results & Performance

### Takeaways & Next Steps


<details><summary>Bibliography</summary>
  <div markdown = "1">
    
## Bibliography

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

inline 
$$ 
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

