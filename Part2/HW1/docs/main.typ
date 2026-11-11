#import "@preview/unequivocal-ams:0.1.2": ams-article, proof, theorem
#import "@preview/typart:0.2.2": donut


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

== kLFS
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
ends are connected. That is, we may visualize the strip as a ring or donut.

#figure(
  scale(50%)[
    #donut(
      (([], rgb("#000")), ([], rgb("#000")), ([], rgb("#DDD")), ([], rgb("#DDD")), ([], rgb("#DDD"))),
      center: $k = 2$ ,
    )
  ],
)


== Experiment 


= Background
#lorem(40)
