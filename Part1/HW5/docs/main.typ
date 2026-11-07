#import "@preview/basic-report:0.5.0": *

#show: it => basic-report(
  heading-color: rgb("#73243d"),
  doc-category: "Homework 5",
  doc-title: "An Application of Takahama's Constrained Optimization Technique",
  author: "Pedro D. Llerenas",
  affiliation: "Centro de Investigación en Matemáticas",
  logo: image("cimat.png", width: 2cm),
  language: "en",
  it,
)

= Introduction

Constrained optimization deals with problems in which we seek the solution to
the problem
$
  min_(x in RR^d) f(x)
$
subject to the following sets
$
  G & = {x in RR^d : g_i (x) = c_i}_(i=1)^n \
  H & ={x in RR^d : h_j (x) <= d_j}_(j = 1)^m
$
One can solve these problems using analytics methods, such as Lagrange
multipliers or KKT conditions. However, one may also take a stochastic
approach to solve the problem. One may observe that the set $G$ has
measure 0 in $RR^d$, thus unreasonable to expect a stochastic process to
land exactly where our equality constraint asks us to.

To bypass this problem, we replace the equality constraints by inequalities,
essentially allowing for an $epsilon$ error around each $c_i$. That is,
$
  G_epsilon = {x in RR^d : norm(g_i (x)) <= epsilon }_(i=1)^n
$
This new set has positive measure, so we can expect a random process to satisfy
this condition with positive probability.

In this document, we perform the analytic calculations of the method.

= Takahama Method in $RR^2$
Let $x in RR^d$ be a feasible point. The iteration depends on the violations.

1. If only 1 restriction $h_i(x) > 0$ is violated: we move in the opposite direction of the gradient of $h_i$. That is,
$
  Delta x = -(h_i (x))/(norm(nabla h_i (x))^2) nabla h_i (x)
$
2. If both restrictions are violated: the step is given by
$
  Delta x = -J^(-1) V,
$
where
$
  J = mat(delim: "[", nabla h_1 (x)^top; nabla h_2 (x)^top)
$
and
$
  V = mat(delim: "[", h_1 (x); h_2 (x)).
$

= Example
Let the objective function be
$
  f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
$
subject to
$
  h_1 (x,y) & = 26 - (x-5)^2 - y^2 & >= 0, \
  h_2 (x,y) & = 20 - 4x - y        & >= 0.
$
We reverse the inequalities to obtain the following restrictions
$
  h_1 (x,y) & = y^2 + (x-5)^2 - 26 & <= 0, \
  h_2 (x,y) & = y + 4x - 20        & <= 0.
$
Now, the violation occurs when $g_i(x,y) > 0$. The gradients of the
restrictions are given by
$
  nabla h_1 (x,y) & = mat(2(x-5); 2y) \
  nabla h_2 (x,y) & = mat(4; 1)
$
== Point 1
We use start at the point $P = (10, 5)$. We evaluate the restrictions at this point,
and obtain the following:
$
  h_1 (10, 5) & = (10 - 5)^2 + 5^2 - 26 & = 24 > 0 \
  h_2 (10, 5) & = 4(10) + 5 - 20        & = 25 > 0
$
Thus, both restriction are violated. The gradients are given by
$
  nabla h_1 (10, 5) & = mat(10; 10) \
  nabla h_2 (10, 5) & = mat(4; 1) \
$
Therefore, the Jacobian evaluates to
$
  abs(J) = mat(delim: "|", 10, 10; 4, 1) = -30,
$
and thus its inverse is
$
  J^(-1) = -1/30 mat(1, -10; -4, 10)
$
so
$
  Delta x = - mat(-1/30, 1/3; 4/30, -1/3) vec(24, 25) = 1/15 vec(-113, 77) approx vec(-7.533, 5.133)
$
Our next evaluation point is then
$
  P_k = P_(k-1) + Delta x = vec(10, 5) + vec(-7.533, 5.133) = vec(2.467, 10.133)
$

== Point 2
Let $P = (5,5)$, then
$
  h_1 (5, 5) & = (5 - 5)^2 + 5^2 - 26 & = -1 <= 0 \
  h_2 (5, 5) & = 4(5) + 5 - 20        &   = 5 > 0
$
Thus, only $h_2$ is violated, so we move in the opposite direction of
its gradient.
$
  nabla h_2 (5, 5) & = mat(4; 1) ==> norm(nabla h_2)&^2 = 17.
$
Therefore,
$
  Delta x = -5/17 vec(4, 1) approx vec(-1.176, -0.294)
$
and the new point is
$
  P_k = P_(k-1) + Delta x = vec(5, 5) + vec(-1.176, -0.294) = vec(3.824, 4.706)
$

== Point 3
Let $P = (2,8)$, then
$
  h_1 (2, 8) & = 47 > 0 \
  h_2 (2, 8) & = -4 <= 0
$
Thus, only $h_1$ is violated, so we move in the opposite direction of
its gradient.
$
  nabla h_1 (2, 8) & = mat(-6; 16) ==> norm(nabla h_1)&^2 = 292.
$
Therefore,
$
  Delta x = -47/292 vec(-6, 16) approx vec(0.966, -2.575)
$
and the new point is
$
  P_k = P_(k-1) + Delta x = vec(2, 8) + vec(0.966, -2.575) = vec(2.966, 5.425)
$

== Point 4
Let $P = (8,8)$, then
$
  h_1 (8, 8) & = 47 > 0 \
  h_2 (8, 8) & = 20 > 0
$
Both restrictions are violated. Thus, the Jacobian is
$
  abs(J) = mat(delim: "|", 6, 16; 4, 1) = -58
$
and thus
$
  J^(-1) = -1/58 mat(1, -16; -4, 6).
$
Therefore,
$
  nabla x = 1/58 vec(-273, -68) approx vec(-4.707, -1.172)
$
and the new point is
$
  P_k = P_(k-1) + Delta x = vec(8, 8) + vec(-4.707, -1.172) = vec(3.293, 6.828).
$
== Graph

#figure(
  image("takahama.png"),
  caption: [Plot of the restrictions and where the points end up after one
  iteration of the Takahama process. We observe that they generally approach
  the feasible lines.],
)


