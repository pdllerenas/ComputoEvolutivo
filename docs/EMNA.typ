#import "@preview/splendid-mdpi:0.1.0"
#import "@preview/lovelace:0.3.1": *


#show: splendid-mdpi.template.with(
  title: [Estimation Multivariate Normal Algorithm with Thresheld Convergence],
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
    month: "August",
    day: 31,
  ),
  keywords: (
    "EMNA",
    "EDA",
    "Genetic Algorithms",
  ),
  venue: image("cimat.png", height: 1.8cm),
  publisher: [],
  paper-type: [],
  venue-link: "www.cimat.mx",
  venue-abbrv: [],
  abstract: [

  ],
  details: none,
)

= Genetic Algorithm (GA) @wikipedia_evolutionary_algorithm
A *genetic algorithm* is a meta heuristic inspired by the process of natural
selection. Genetic algorithms are used to generate (usually) non-optimal but
high-quality solutions to optimization and search problems via biologically
inspired operators such as _selection_, _crossover_, and _mutation_.

= Estimation Distribution Algorithm (EDA)
An *estimation distribution algorithm* is a meta heuristic derived from a GA,
where the operations of crossover and crossover are replaced by a probability
distribution function (PDF), which is determined by the fittest individuals in
the population. This model is used to replace a portion (0% to 100%) of the
population. And thus, and EDA converts an optimization problem into finding a
fitting PDF for our objective function, in which we except the PDF to have
a higher density near or at local minima.

= Estimation Multivariate Normal Algorithm (EMNA)
An *estimation multivariate normal algorithm* is an EDA in which the proposed
distribution used to generate the novel populations is a multivariate normal (MVN).
We use the fittest individuals mean, variance and covariance to determine the
parameters of the MVN. The size of the sample used is given by the _selection
coefficient (sc)_, which is usually set to 0.30.

#figure(
  kind: "algorithm",
  supplement: [Algorithm],

  pseudocode-list(booktabs: true, numbered-title: [EMNA])[
    + $X :=$ Initial Population
    + *while* evals < maxEvals
      + $X_s :=$ Selection($X$, sc)
      + $[mu, Sigma] :=$ ParameterEstimation($X_s$)
      + $X :=$ SampleMVN($mu, Sigma$)
    + *end*
    + *return* bestSolution
  ],
)


= EMNA with Thresheld Convergence @tamayo2015emna
An *EMNA with Thresheld Convergence* proposes using a minimum step size at each generation
for the newer solutions. Since there is no parent-child relationship in an EMNA,
the step is applied on the parameter space. Specifically, we restrict the norm
of the covariance matrix generated from the elite individuals. We define the _umbral_ as
$
  "umbral" = alpha times "diagonal" times (("maxEvals" - "evals")/"maxEvals")^gamma.
$
We use this number to perform the following action on the covariance matrix ($Sigma$):
$
  Sigma := "umbral" * Sigma / norm(Sigma).
$
This prevents the EMNA from prematurely converging.

= Results
We used 3 test functions @surjanovic_griewank, @surjanovic_ackley, @surjanovic_sphere on the algorithm, for $d = 2, 5$:
$
    "Sphere: " & f(x) = sum_(i = 1)^d x_i^2 \
  "Griewank: " & f(x) = sum_(i=1)^d (x_i^2)/4000 - product_(i=1)^d cos(x_i/sqrt(i)) + 1 \
    "Ackley: " & f(x) = -a exp(- b sqrt(1/d sum_(i=1)^d x_i^2)) - exp(1/d sum_(i=1)^d cos(c x_i)) + a +exp(1) \
$
We used a population of $n = 1000$, with a maximum of $100000 times d$
functions evaluations (one more zero than in the original paper).

#figure(
  table(
    stroke: none,
    columns: 6,
    table.hline(),
    table.header([*Function*], [*Best*], [*Worst*], [*Mean*], [*Median*], [*StdDev*]),
    table.hline(),
    [Ackley],
    [$9.97 times 10^(-3)$],
    [$1.16 times 10^(-1)$ ],
    [$4.63 times 10^(-2)$],
    [$3.92 times 10^(-2)$ ],
    [$2.96 times 10^(-2)$ ],
    [Griewank],
    [$3.44 times 10^(-3)$],
    [$1.50 times 10^(-2)$ ],
    [$8.50 times 10^(-3)$],
    [$8.59 times 10^(-3)$ ],
    [$3.63 times 10^(-3)$ ],
    [Sphere],
    [$7.46 times 10^(-6)$],
    [$9.79 times 10^(-5)$ ],
    [$6.39 times 10^(-5)$],
    [$7.86 times 10^(-5)$ ],
    [$3.16 times 10^(-5)$ ],
    table.hline(),
  ),
  caption: [Statistics of EMNA for 10 randomly seeded runs on $d = 2$ ],
)


#figure(
  table(
    stroke: none,
    columns: 6,
    table.hline(),
    table.header([*Function*], [*Best*], [*Worst*], [*Mean*], [*Median*], [*StdDev*]),
    table.hline(),
    [Ackley],
    [$5.05 times 10^(-1)$],
    [$1.79 times 10^(0)$],
    [$1.13 times 10^(0)$],
    [$1.27 times 10^(0)$],
    [$4.35 times 10^(-1)$],
    [Griewank],
    [$1.64 times 10^(-1)$],
    [$4.38 times 10^(-1)$],
    [$2.52 times 10^(-1)$],
    [$2.52 times 10^(-1)$],
    [$7.24 times 10^(-2)$],
    [Sphere],
    [$2.86 times 10^(-4)$],
    [$2.12 times 10^(-3)$],
    [$1.37 times 10^(-3)$],
    [$1.42 times 10^(-3)$],
    [$5.03 times 10^(-4)$],
    table.hline(),
  ),
  caption: [Statistics of EMNA for 10 randomly seeded runs on $d = 5$],
)

= Conclusions
We observed that EMNA struggled to find the true minimum of the Ackley function
for both test cases. The other two functions we closer to converging. This
makes sense, as the Ackley function is nearly flat outside the origin, with
many local minima, while the origin is a steep drop-off to the true global
minimum. The Griewank function has the form of a paraboloid in 2D, where
locally, the curve looks like a cosine function in 2D. This leads to many local
minima, while the true global minimum is in the origin. Compared to Ackley,
this functions local minima converge towards the global minimum, although it is
possible to get stuck in a local minimum. Finally, the sphere function is quite
simple, allowing for a fast convergence to the global minimum.

#bibliography("refs.bib")
