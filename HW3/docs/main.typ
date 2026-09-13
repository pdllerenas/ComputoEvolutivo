#import "@preview/splendid-mdpi:0.1.0"
#import "@preview/lovelace:0.3.1": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, edge, node

// #set page(fill: rgb("#1e1e1e"))
// #set text(fill: rgb("#e0e0e0"))

#show: splendid-mdpi.template.with(
  title: [Homework 3: KL Divergence and Maximum Likelihood Estimation of the Graphic Model],
  authors: (
    (
      name: "Pedro Llerenas",
      department: "MCC",
      institution: "CIMAT",
      city: "Guanajuato",
      country: "Mexico",
      mail: "pedro.llerenas@cimat.mx",
    ),
  ),
  date: (
    year: 2026,
    month: "September",
    day: 22,
  ),
  keywords: (
    "Kullback-Leibler Divergence",
    "Mutual Information",
    "Normal Distribution",
    "Maximum Likelihood",
    "Graphic Model",
  ),
  venue: image("cimat.png", height: 1.8cm),
  publisher: [],
  paper-type: [],
  venue-link: "www.cimat.mx",
  venue-abbrv: [],
  abstract: [
    We explore three models for two databases. One model considers independent
    random variables, another one is a chain of dependencies, and the last one
    is a dependency tree. We calculate their maximum likelihood and their
    pair-wise Kullback-Leibler divergence, where we observe...
  ],
  details: none,
)


#let dkl(p, q) = $D_("KL")(#p parallel #q)$

= Introduction
Being able to predict an outcome of a certain event will remain to be the most
powerful tool humanity possesses. To be able to do this, we attempt to describe
these events with mathematical models. Although these models cannot accurately
predict a future events with 100% certainty, being able to deduce that an
outcome is _likely_ to occur, allows judging the situations an easier task.

Once an event is has been sampled over many iterations, we can observe
patterns in the data they generate. We use known probability distributions
that may follow these patterns. To be able to justify why a certain distribution
_fits_ the observed data, we use maximum likelihood estimates, which tell us
how likely it is to observe what was observed, if we assume a certain distribution
is the source of the data. Maximizing this likelihood then leads to properly adjusted
distributions.

In this document, we use Mutual-Information-Maximizing Input Clustering

= Preliminaries
Before we present the models, we provide some definitions.

== Entropy
For a discrete random variable $X$, we denote by $H(X)$ the _entropy_ of $X$. We define it as
$
  H(X) = -sum_(x in X) P(x) log P(x).
$
This value quantifies the average level of uncertainty or information
associated with the variable's potential states or possible outcomes.

== Mutual Information
For a pair of random variables $X,Y$, the mutual information $I(X;Y)$ measures
how much information about $X$ is gained from knowing $Y$. It is defined as
$
  I(X;Y) = sum_(x in X) sum_(y in Y) P(x, y) log P(x mid(|) y)/P(x)
$

== Maximum Likelihood Estimation
This method estimates the parameters of an assumed probability distribution,
given some observed data. We model a set of observations as a random sample
from an unknown joint probability distribution which is expressed in terms of a
set of parameters. We write the parameters governing the joint distribution as
a vector $theta = [theta_1, ..., theta_k]^top$ so that this distribution falls
within a parametric family ${f(dot;theta) | theta in Theta}$, where $Theta$ is
called the parameter space, a finite-dimensional subset of Euclidean space.
Evaluating the joint density at the observed data sample $bold(y) = (y_1, ...,
  y_n)$ gives a real-valued function,
$
  cal(L)_n (theta) = cal(L)_n (theta;bold(y)) = f_n (bold(y); theta)
$
which is called the likelihood function. For independent random variables $f_n
(bold(y);theta)$ will be the product of univariate density functions:
$
  f_n (bold(y);theta) = product_(k=1)^n f_k^("univar") (y_k; theta).
$
Consequently, the log-likelihood is given by
$
  log(f_n (bold(y);theta)) = sum_(k=1)^n log f_k^("univar") (y_k; theta).
$

= Models
== Independent Variables
We first explore the simplest model possible: independent random variables.
This produces a joint probability distribution
$
  P(X) = product_(1 <= i <= n) P(X_i)
$
== MIMIC
This algorithm was first proposed in @DeBonet1996MIMICFO, it takes a
combinatorial approach to the best model. That is, it produces a first-order
dependency chain (a simple Bayesian network) to represent the joint
probability. It finds the optimal ordering $X_(pi(1)), X_(pi(2)), ...,
X_(pi(n))$ by minimizing the Kullback-Leibler divergence between the true
distribution of the elite solutions and the estimated chain distribution.
Mathematically, this reduces to finding the chain that maximizes the mutual
information between adjacent variables.


#figure(
  kind: "algorithm",
  supplement: [Algorithm],

  pseudocode-list(booktabs: true, numbered-title: [MIMIC])[
    + $X :=$ Random initial population of size $N$
    + Evaluate fitness $f(x)$ for all $x in X$
    + *while* not converged *do*
      + $theta_t :=$ Fitness threshold for top $theta$-percentile
      + $E := {x in X mid(|) f(x) >= theta_t}$ (Elite subset)
      + Compute marginal $H(X_i)$ and conditional $H(X_i | X_j)$ entropies for all pairs in $E$
      + Find permutation $pi$ to minimize total conditional entropy: $H(X_(pi(1))) + sum_(i=1)^(n-1) H(X_(pi(i+1)) | X_(pi(i)))$
      + $X :=$ Sample $N$ new candidates from this dependency chain
      + Evaluate fitness $f(x)$ for all new $x in X$
    + *end*
    + *return* Best candidate $x^ast in X$
  ],
)

In this document, we will focus on $theta = 100%$. That is, without any elitism.
The algorithm then simplifies to the following.

#figure(
  kind: "algorithm",
  supplement: [Algorithm],

  pseudocode-list(booktabs: true, numbered-title: [Simplified MIMIC])[
    + $X :=$ Random initial population of size $N$
    + Compute marginal $H(X_i)$ and conditional $H(X_i | X_j)$ entropies for all pairs.
    + Find permutation $pi$ to minimize total conditional entropy: $H(X_(pi(1))) + sum_(i=1)^(n-1) H(X_(pi(i+1)) | X_(pi(i)))$
    + *return* $pi$
  ],
)

#figure(
  diagram(
    node-stroke: .1em,
    node-fill: gradient.radial(white.lighten(80%), gray.lighten(40%), center: (20%, 20%), radius: 20%),
    spacing: 4em,
    node((0, 0), $x_1$, radius: 1.5em),
    edge("-|>"),
    node((1, 0), $x_3$, radius: 1.5em),
    edge("-|>"),
    node((2, 0), $x_2$, radius: 1.5em),
    edge("-|>"),
    node((3, 0), $x_4$, radius: 1.5em),
  ),
  caption: [Example of a first-order dependency chain.],
)
== Dependency Tree
This model is similar to MIMIC in their goal of maximizing mutual information
(or minimize entropy). However, dependency tree does not restrict the number of
children a node can have, producing a directed acyclic graph, where each
variable has at most one parent.

The following algorithm is also a simplified version of dependency tree, where
we do not consider elitism. It was first presented in @ChowLiu.
#figure(
  kind: "algorithm",
  supplement: [Algorithm],

  pseudocode-list(booktabs: true, numbered-title: [Dependency Tree])[
    + $X :=$ Random initial population of size $N$
    + Compute Mutual Information $I(X_i;X_j)$ for all pairs of variables
    + Let $G = (V, E)$ be a complete graph where weight $W(i,j) = I(X_i; X_j)$
    + $T :=$ Maximum Spanning Tree of $G$ (via Kruskal)
    + Select an arbitrary node $R in V$ as the root
    + Direct all edges in $T$ outward from $R$ to form a dependency tree
    + Estimate marginal $P(R)$ and conditionals $P(X_i mid(|) X_("parent"(i)))$
    + *return* pointer to $R$ (root of tree)
  ],
)

= Experiment Setup
== Experiment 1
We will use be using 6 binary random variables for which we have the following empirical data
#figure(
  table(
    columns: 4,
    align: center + horizon,
    fill: (x, y) => if calc.odd(y) { luma(200) } else { white },
    stroke: none,
    table.hline(),
    table.header([*Variable* ($X_i$)], [*Parent* ($X_(i-1)$)], [$P(X_i=1 | X_(i-1)=0)$], [$P(X_i=1 | X_(i-1)=1)$]),
    table.hline(stroke: .5pt),
    $x_1$, [None], [--], [0.20],
    $x_2$, $x_1$, [0.92], [0.31],
    $x_3$, $x_2$, [0.42], [0.15],
    $x_4$, $x_3$, [0.28], [0.75],
    $x_5$, $x_4$, [0.67], [0.43],
    $x_6$, $x_5$, [0.35], [0.82],
    table.hline(),
  ),
  caption: [Conditional probabilites for the first database DB1.],
)

#figure(
  diagram(
    node-stroke: .1em,
    node-fill: gradient.radial(white.lighten(80%), gray.lighten(40%), center: (20%, 20%), radius: 20%),
    spacing: 4em,
    node((0, 0), $x_1$, radius: 1.5em),
    edge("-|>"),
    node((1, 0), $x_2$, radius: 1.5em),
    edge("-|>"),
    node((2, 0), $x_3$, radius: 1.5em),
    edge("-|>"),
    node((3, 0), $x_4$, radius: 1.5em),
    edge("-|>"),
    node((4, 0), $x_5$, radius: 1.5em),
    edge("-|>"),
    node((5, 0), $x_6$, radius: 1.5em),
  ),
  caption: [Dependency chain from DB1. Its joint probability is given by $P(X) = P(X_1)product_(i = 2)^6 P(X_i|P_(i-1))$ ],
)

This data resembles a simple dependency chain, for which we expect MIMIC to
perform well. Dependency tree should also perform well.

== Experiment 2
The next experiment also contains 6 binary random variables.
#figure(
 table(
    columns: 4,
    align: center + horizon,
    fill: (x, y) => if calc.odd(y) { luma(200) } else { white },
    stroke: none,
    table.hline(),
    table.header(
      [*Variable* ($X_i$)], [*Parent* ($X_"parent"$)], [$P(X_i=1 | X_"parent"=0)$], [$P(X_i=1 | X_"parent"=1)$]
    ),
    table.hline(stroke: .5pt),
    $x_1$, [None], [--], [0.20],
    $x_2$, $x_1$, [0.92], [0.31],
    $x_4$, $x_1$, [0.42], [0.15],
    $x_5$, $x_2$, [0.28], [0.75],
    $x_6$, $x_2$, [0.67], [0.43],
    $x_3$, $x_6$, [0.35], [0.82],
    table.hline(),
  ),
  caption: [Conditional probabilites for the second database DB2.],
)
In this case, $x_1$ contains 2 children, and so does $x_2$. We expect
dependency tree to perform well, but MIMIC should not be able to reproduce the
same dependency list (by design).

#figure(
  diagram(
    node-stroke: .1em,
    node-fill: gradient.radial(white.lighten(80%), gray.lighten(40%), center: (20%, 20%), radius: 20%),
    spacing: 4em,
    node((0, 0), $x_1$, radius: 1.5em),

    edge((0, 0), (1, 1), "-|>"),
    edge((0, 0), (-1, 1), "-|>"),

    node((-1, 1), $x_2$, radius: 1.5em),
    node((1, 1), $x_4$, radius: 1.5em),

    edge((-1, 1), (-1.5, 2), "-|>"),
    edge((-1, 1), (-0.5, 2), "-|>"),

    node((-0.5, 2), $x_5$, radius: 1.5em),
    node((-1.5, 2), $x_6$, radius: 1.5em),

    edge((-1.5, 2), (-1.75, 3), "-|>"),

    node((-1.75, 3), $x_3$, radius: 1.5em),
  ),
  caption: [Dependency tree generated by DB2. Its joint probability is given by $P(X) = P(X_1)P(X_2|X_1)P(X_4|X_1)P(X_6|X_2)P(X_5|X_2)P(X_3|X_6)$ ],
)

#bibliography("refs.bib")
