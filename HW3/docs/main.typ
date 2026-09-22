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
  I(X;Y) = sum_(x in X) sum_(y in Y) P(x, y) log P(x mid(|) y)/P(x) = H(X) - H(X|Y) = H(Y) - H(Y|X)
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
  caption: [Conditional probabilities for the second database DB2.],
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

= Results

== Independent Variables (Model 1)

=== Maximum Likelihood Estimation
For this model, since we do not have dependencies, we only need to observe the marginal probabilities.
For an execution with a database size of 30,000 we obtained the following probabilities, which also
correspond to the maximum likelihood estimation.
#grid(
  columns: 2,
  figure(
    table(
      columns: 2,

      stroke: none,
      align: center + horizon,
      fill: (x, y) => if calc.odd(y) { luma(200) } else { white },
      table.hline(),
      table.header([Bit], [Probability]),
      table.hline(),
      $P(X_1 = 1)$, $0.20198$,
      $P(X_2 = 1)$, $0.79633$,
      $P(X_3 = 1)$, $0.20297$,
      $P(X_4 = 1)$, $0.37415$,
      $P(X_5 = 1)$, $0.58213$,
      $P(X_6 = 1)$, $0.62222$,
      table.hline(),
    ),
    caption: [Values for Database 1.],
  ),

  figure(
    table(
      columns: 2,

      align: center + horizon,
      fill: (x, y) => if calc.odd(y) { luma(200) } else { white },
      stroke: none,
      table.hline(),
      table.header([Bit], [Probability]),
      table.hline(),
      $P(X_1 = 1)$, $0.2023$,
      $P(X_2 = 1)$, $0.79675$,
      $P(X_3 = 1)$, $0.22539$,
      $P(X_4 = 1)$, $0.37619$,
      $P(X_5 = 1)$, $0.48014$,
      $P(X_6 = 1)$, $0.72331$,
      table.hline(),
    ),
    caption: [Values for Database 2.],
  ),
)


Despite having identical probabilities a priori, the values are different, as
we used the same number generator to first generate DB1, then DB2, without
resetting it. This is the case throughout the whole program; we only use 1
instance of `std::mt19937`.

=== Kullback-Leibler Divergence
We obtained $D_("KL")(M_1 || "DB1") = 0.587202$, and
$D_("KL")(M_1 || "DB2") = 0.587202$.


== MIMIC (Model 2)
=== Resulting Chain
We produced the chain
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
  caption: [Chain generated by MIMIC for DB1.],
)
Which corresponds to the original model used to generate the database. Note
that sometimes (rarely, maybe $p < 0.05$), we obtain slightly different chains. This is due
to the randomness of the generation, which can cause some weakly dependent
nodes to lose their original pairing.

// 1 -> 2 -> 6 -> 3 -> 4 -> 5
For database 2, the following chain was generated.
#figure(
  diagram(
    node-stroke: .1em,
    node-fill: gradient.radial(white.lighten(80%), gray.lighten(40%), center: (20%, 20%), radius: 20%),
    spacing: 4em,
    node((0, 0), $x_1$, radius: 1.5em),
    edge("-|>"),
    node((1, 0), $x_2$, radius: 1.5em),
    edge("-|>"),
    node((2, 0), $x_6$, radius: 1.5em),
    edge("-|>"),
    node((3, 0), $x_3$, radius: 1.5em),
    edge("-|>"),
    node((4, 0), $x_4$, radius: 1.5em),
    edge("-|>"),
    node((5, 0), $x_5$, radius: 1.5em),
  ),
  caption: [Chain generated by MIMIC for DB2.],
)
It was expected to not even be close to the true nature of the database, as it was
originally a tree, so we cannot expect a simple chain to approximate it well.

=== Maximum Likelihood Estimation
The produced joint probability values are presented in the following tables.
#grid(
  columns: 2,
  figure(
    table(
      columns: 2,
      stroke: none,
      align: center + horizon,
      fill: (x, y) => if calc.odd(y) { luma(200) } else { white },
      table.hline(),
      table.header([Bit], [Probability]),
      table.hline(),
      $P(x_1 = 1)$, $0.20198$,
      $P(x_2 = 1 | x_1 = 0)$, $0.92123$,
      $P(x_2 = 1 | x_1 = 1)$, $0.302852$,
      $P(x_3 = 1 | x_2 = 0)$, $0.419158$,
      $P(x_3 = 1 | x_2 = 1)$, $0.147677$,
      $P(x_4 = 1 | x_3 = 0)$, $0.277769$,
      $P(x_4 = 1 | x_3 = 1)$, $0.752624$,
      $P(x_5 = 1 | x_4 = 0)$, $0.672158$,
      $P(x_5 = 1 | x_4 = 1)$, $0.431538$,
      $P(x_6 = 1 | x_5 = 0)$, $0.348194$,
      $P(x_6 = 1 | x_5 = 1)$, $0.818924$,
      table.hline(),
    ),
    caption: [Probability tables by MIMIC for DB1.],
  ),

  figure(
    table(
      columns: 2,
      stroke: none,
      align: center + horizon,
      fill: (x, y) => if calc.odd(y) { luma(200) } else { white },
      table.hline(),
      table.header([Bit], [Probability]),
      table.hline(),
      $P(x_1 = 1)$, $0.2023$,
      $P(x_2 = 1 | x_1 = 0)$, $0.920133$,
      $P(x_2 = 1 | x_1 = 1)$, $0.310232$,
      $P(x_6 = 1 | x_2 = 0)$, $0.3477$,
      $P(x_6 = 1 | x_2 = 1)$, $0.819128$,
      $P(x_3 = 1 | x_6 = 0)$, $0.419567$,
      $P(x_3 = 1 | x_6 = 1)$, $0.151111$,
      $P(x_4 = 1 | x_3 = 0)$, $0.36643$,
      $P(x_4 = 1 | x_3 = 1)$, $0.409734$,
      $P(x_5 = 1 | x_4 = 0)$, $0.461791$,
      $P(x_5 = 1 | x_4 = 1)$, $0.510566$,
      table.hline(),
    ),
    caption: [Probability tables by MIMIC for DB2.],
  ),
)

=== KL-Divergence
We obtained $D_("KL")(M_2 || "DB1") = 0.000453911$, and
$D_("KL")(M_1 || "DB2") = 0.587202$.

== Chow-Liu Tree
=== Tree Generated
#figure(
  diagram(
    node-stroke: .1em,
    node-fill: gradient.radial(white.lighten(80%), gray.lighten(40%), center: (20%, 20%), radius: 20%),
    spacing: 4em,

    // Root node
    node((0, 0), $x_5$, radius: 1.5em),

    // Edges from root (x5 -> x4, x5 -> x6)
    edge((0, 0), (-1, 1), "-|>"),
    edge((0, 0), (1, 1), "-|>"),

    // Level 1 nodes
    node((-1, 1), $x_4$, radius: 1.5em),
    node((1, 1), $x_6$, radius: 1.5em),

    // Edge x4 -> x3
    edge((-1, 1), (-1, 2), "-|>"),

    // Level 2 node
    node((-1, 2), $x_3$, radius: 1.5em),

    // Edge x3 -> x2
    edge((-1, 2), (-1, 3), "-|>"),

    // Level 3 node
    node((-1, 3), $x_2$, radius: 1.5em),

    // Edge x2 -> x1
    edge((-1, 3), (-1, 4), "-|>"),

    // Level 4 node
    node((-1, 4), $x_1$, radius: 1.5em),
  ),
  caption: [Dependency tree generated by the Chow-Liu Algorithm for DB1.],
)
We observe that although it is not exactly the chain the database came from, this is expected,
as taking an arbitrary root node leads to this. We can observe that if we instead take node 1
as the root, it leads to the original chain (by reversing the edges).

#figure(
  diagram(
    node-stroke: .1em,
    node-fill: gradient.radial(white.lighten(80%), gray.lighten(40%), center: (20%, 20%), radius: 20%),
    spacing: 4em,

    // Root node
    node((0, 0), $x_1$, radius: 1.5em),

    // Edges from root (x1 -> x2, x1 -> x4)
    edge((0, 0), (-1, 1), "-|>"),
    edge((0, 0), (1, 1), "-|>"),

    // Level 1 nodes
    node((-1, 1), $x_2$, radius: 1.5em),
    node((1, 1), $x_4$, radius: 1.5em),

    // Edges from x2 (x2 -> x5, x2 -> x6)
    edge((-1, 1), (-2, 2), "-|>"),
    edge((-1, 1), (0, 2), "-|>"),

    // Level 2 nodes
    node((-2, 2), $x_5$, radius: 1.5em),
    node((0, 2), $x_6$, radius: 1.5em),

    // Edge from x6 -> x3
    edge((0, 2), (0, 3), "-|>"),

    // Level 3 node
    node((0, 3), $x_3$, radius: 1.5em),
  ),
  caption: [Dependency tree generated by the Chow-Liu Algorithm for DB2.],
)
This tree is verbatim to the original tree. Since it is random, it is also
possible to obtain other permutations of the tree, where the root may not be
node 1. However, the joint probability stays the same.

=== Maximum Likelihood Estimate
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
    $x_5$, [None], [--], [0.58213],
    $x_4$, $x_5$, [0.508986], [0.277361],
    $x_3$, $x_4$, [0.0802269], [0.408285],
    $x_2$, $x_3$, [0.851574], [0.579396],
    $x_1$, $x_2$, [0.691363], [0.0768149],
    $x_6$, $x_5$, [0.348194], [0.818924],
    table.hline(),
  ),
  caption: [Conditional probabilities for the second database DB1 generated by the Chow-Liu tree.],
)

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
    $x_1$, [None], [--], [0.2023],
    $x_2$, $x_1$, [0.920133], [0.310232],
    $x_5$, $x_2$, [0.672226], [0.431139],
    $x_6$, $x_2$, [0.3477], [0.819128],
    $x_3$, $x_6$, [0.419567], [0.151111],
    $x_4$, $x_1$, [0.280343], [0.754128],
    table.hline(),
  ),
  caption: [Conditional probabilities for the second database DB2 generated by the Chow-Liu tree.],
)

=== KL-Divergence
We obtained $D_("KL")(M_2 || "DB1") = 0.000453911$, and
$D_("KL")(M_1 || "DB2") = 0.000399334$.

= Conclusions
We may conclude that Chow-Liu tree produces the best results for both databases in general. In particular, it does not suffer
the 1-dimensionality problem of the MIMIC model, which only allows for 1 child node. The independent nodes model is extremely simple,
but it does not capture the true nature of the databases, since it completely ignores the joint probabilities.

= Ideas
One possible way of allowing the MIMIC model to have greater versatility is to
choose 1 root node, and allow only the root node to have various children. This
keeps an almost linear model, but may allow leeway for some generalizations.

= Appendix

== Compilation Instructions
Requirements:
- C++20
- CMake
- Make

First, `cd` into the directory. Then,

```sh
mkdir build && cd build
```
Then, we use CMake to generate the make files.

```sh
cmake ..
```
Finally, we use `make` to compile and obtain the necessary binary files.
```sh
make
```

== Execution Instructions
All the results are obtained from a single execution, all printed to `cout`. Once 
all compilations are done, we do (while in the `build` directory)
```sh
./src/main
```
This will generate results for all 3 models, for both databases.


#bibliography("refs.bib")
