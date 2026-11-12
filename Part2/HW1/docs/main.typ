#import "@preview/unequivocal-ams:0.1.2": ams-article, proof, theorem
#import "@preview/typart:0.2.2": donut

#show link: underline
#show link: set text(fill: blue)

#show: ams-article.with(
  title: [Solving the Frequency Assignment Problem with Cellular Algorithms],
  authors: (
    (
      name: "Pedro D. Llerenas",
      department: [Maestría en Ciencias de la Computación],
      organization: [Centro de Investigación en Matemáticas],
      location: [Guanajuato, Guanajuato, MX],
      email: "pedro.llerenas@cimat.mx",
      url: "",
    ),
  ),
  abstract: [This paper aims to solve the Frequency Assignment Problem (FAP)
    through a cellular algorithm, namely, k-Fixed Line Sweep (kFLS)
    . We observe the effect of $k$ over its whole range, and find that // TODO: explain
  ],
  bibliography: bibliography("refs.bib"),
)

= Introduction
The Frequency Assigment Problem (FAP) deals with the issue of the growing
number of wireless communication systems deployed around the globe, with only a
limited radio frequency spectrum. This affects planning models for permanent
spectrum allocation, licensing, regulation and network design. Furthermore, the
algorithms for dynamically assigning frequencies to users within an established
network @murphey1999. This paper takes a stochastic approach to the problem in
form of a cellular algorithm. Specifically, we will be using kFLS
@cornejo2021first.

= k-LFS
The scheme consists of the usual Fixed Line Sweep (FLS) with a
slight modification, as follows:

1. Synchronously explore $k$ adjacent cells instead of exploring cells one by one.
2. Asynchronously explore the cells by shifting the $k$ adjacent cells by one to the right.

The special cases of $k = 1$ and $k = n$, where $n$ is the population size,
correspond to FLS and synchronous exploration of the whole population,
respectively.

#figure(
  grid(
    columns: 5,

    // By a single value
    align: center,
    // By a single but more complicated value
    inset: (x: 8pt, y: 8pt),
    // By an array of values (cycling)
    fill: (x, y) => if calc.rem(6 - x + y, 5) < 2 { black } else { rgb("#DDD") },
    // By a function that returns a value
    stroke: 0.5pt,

    ..range(5 * 5).map(n => none)
  ),
  caption: [kLFS for $k = 2$],
)

We make the observation that the strip has the topology of a torus, in the sense that the
ends are connected. That is, we may visualize the strip as a ring or donut, as depicted in @donut.

#figure(
  box(
    width: 100%,
    height: 35%,
    clip: true,
    scale(50%)[
      #donut(
        (([], rgb("#000")), ([], rgb("#000")), ([], rgb("#DDD")), ([], rgb("#DDD")), ([], rgb("#DDD"))),
        center: $k = 2$,
      )
    ],
  ),
  caption: [Visualization of topology of East-West update policy, for a single update block.],
)<donut>


= Frequency Assignment Problem
As mentioned in @murphey1999, a frequency assignment problem (FAP) models the
task of assigning radio spectrum to a set of transmitters. That is, if $V$  is
the set of all transmitters with $v in V$, $f(v)$ denotes the continuum of
frequencies assigned to $v$ by the assignment rule $f$. A feasible frequency
assignment must additionaly satisfy certain interference and bandwidth
constraints, $I$ and $B$ respectively. An optimal frequency assignment globally
minimizes a cost function that depends upon a particular objective of the
problem. Two possible objectives and thus two basic problems formulations are
widely used:

- minimize a bandwidth function subject to an acceptable level of interference
- minimize the system interference given a fixed frequency allotment

where the bandwidth function is either the span of frequencies (difference
between maximum and minimum frequencies) assigned or the order of channels
assigned. The system interference function is less concise but can be defined
as the number of resourced denied service, or the bit error rate for the
combined network. This problem is known to be NP-hard, meaning that a heuristic approach
is justified. In this paper, we apply a constant penalization as follows:

$
  P(u,v) = cases(
    10 & " if " abs(f_u - f_v) = 1,
    100 & " if " f_u = f_v,
    0 & " otherwise"
  )
$
This an extremely simplified penalization (fitness) function that does not
consider the distance between the nodes. We plan on addressing these lacking
features in a future paper.

== Implementation

We mapped the graph topology into a flat array representation as follows:

- Frequency Domain Boundaries: Valid frequencies are bounded within an integer
  range $[f_min, f_max]$.
- Edge-List Storage: Interference constraints are stored as an unweighted edge
  vector of `pair<int, int>`.

Note that they are left unweighted because we will not address the distance
aspect of the problem in this paper.

= Experiment
We will do experiments on graphs with unweighted edges. That is, no notion of
distance will be present. This simplifies the problem.

To swiftly obtain results, we used MPI to distribute the workload of various indices of
$k$ (for which we execute 10 runs) to 10 distinct processors. We also use
OpenMP to distribute the 10 runs into various threads. We made use of
_CIMAT's Laboratorio de Supercomputo_.

== Instances
The standard instances to use for this problem are from the
#link("https://fap.zib.de/problems/")[FAP Website]. However, these
instances also include additional problem constraints in their
initialization files. Therefore, we opted to instead use instances from
the #link("https://mat.tepper.cmu.edu/COLOR/instances.html")[graph coloring problem].
We will be using 3 specific files:
- `fpsol2.i.1.col`,
- `fpsol2.i.2.col`,
- `fpsol2.i.3.col`,

The structure of the files is as follows. The first character of the file
determines what the following line of text represents.
#figure(
  ```
  c - COMMENTS
  p - HEADER (dimensions of problem)
  e - EDGES
  ```,
)

For example, `fpsol2.i.2.col` has the line `p edge 451 8691`, which means there
are 451 nodes and a total of 8691 edges.

== Results
On all three instances, we intialize the individuals at with a random frequency
in the integer range ${0, ... ,|V|}$ uniformly. It is possible to find a better initial
estimate to accelerate convergence @Brooks_1941. Their idea involved counting
the number of edges each node has, and obtain its maximum. This gives a tighter upper
bound on the problem. We do no implement this approach.

=== Instance 1
The first instance `fpsol2.i.1.col` has 496 nodes and 11654 edges. According to the
file source page, the optimal color count is 65. That is, 65 distinct colors are needed to
solve the associated Graph Coloring Problem. In our case, we only aim to minimize the
penalty function.

#figure(
  image("images/e1_kfls_combined_plot.png", width: 80%),
  caption: [Solutions found for instance 1. No visible pattern is seen on $k$.],
)

=== Instance 2
The first instance `fpsol2.i.2.col` has 451 nodes and 8691 edges. According to the
file source page, the optimal color count is 30.

#figure(
  image("images/e2_kfls_combined_plot.png", width: 80%),
  caption: [Solutions found for instance 2. No visible pattern is seen on $k$.],
)

=== Instance 3
The first instance `fpsol2.i.3.col` has 425 nodes and 8688 edges. According to the
file source page, the optimal color count is 30.

#figure(
  image("images/e3_kfls_combined_plot.png", width: 80%),
  caption: [Solutions found for instance 3. No visible pattern is seen on $k$.],
)
None of the instances show any dependence on $k$. Moreover, the solutions found
for each $k$ have a high variance, essentially being separated in two groups.
They either stay far away from the solution, or stay within a 600 penalty range.
This could be due to the omission of the distance parameter as a whole, which
might make the algorithm too stochastic; not entirely able to make good decisions.

= Conclusions


#lorem(40)
