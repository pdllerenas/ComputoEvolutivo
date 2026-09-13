
#import "@preview/splendid-mdpi:0.1.0"
#import "@preview/lovelace:0.3.1": *

// #set page(fill: rgb("#1e1e1e"))
// #set text(fill: rgb("#e0e0e0"))

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
    day: 7,
  ),
  keywords: (
    "Kullback-Leibler Divergence",
    "Mutual Information",
    "Normal Distribution",
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


#let dkl(p, q) = $D_("KL")(#p parallel #q)$

= Kullback-Leibler Divergence
The *Kullback-Leibler (KL) divergence*, denoted by $dkl(P, Q)$, measures
how much an approximating probability distribution $Q$ is different from a true
probability distribution $P$. Mathematically, it is defined as
$
  dkl(P, Q) = sum_(x in X) P(x) log(P(x)/Q(x)).
$
One may also extend the definition to continuous random variables as
$
  dkl(P, Q) = integral_RR p(x) log(p(x)/q(x)) dif x
$

An interpretation of the KL divergence of $P$ from $Q$ is the expected excess
surprisal from using the approximation $Q$  instead of $P$ when the actual is
$P$. It is important to note that this is not a metric, as in general, the
symmetric property $dkl(P, Q) = dkl(Q, P)$ does not hold (neither does the triangle
inequality). It does satisfy the properties $dkl(P, P) = 0$ and $dkl(P, Q) >= 0$.

= KL-Divergence of Normal Distributions
A particularly interesting application of the KL divergence is to two normal distributions,
which may pop-up in almost any real-life phenomenon by the central limit theorem. This can help
us gauge how much information is lost when attempting to approximate one normal distribution
with another.

Let $P = cal(N)(mu_1, sigma_1^2)$ and $Q = cal(N)(mu_2, sigma_2^2)$. We first note that
$
  dkl(P, Q) & = integral_RR p(x) log(p(x)/q(x)) dif x \
            & = integral_RR p(x) log p(x) dif x - integral_RR p(x) log q(x) dif x\.
$
The first term gives

$
  integral_RR p(x) log p(x) dif x &= integral_RR (- 1/2 log(2 pi sigma_1^2) -(x-mu_1)^2/(2sigma_1^2)) p(x) dif x\
  &= - 1/2 log(2 pi sigma_1^2) underbrace(integral_RR p(x) dif x, 1) - 1/(2sigma_1^2) underbrace(integral_RR (x-mu_1)^2p(x), sigma_1^2) dif x\
  &= - 1/2 log(2 pi sigma_1^2) - 1/2.
$
For the second term, we have
$
  integral_RR p(x) log q(x) dif x &= integral_RR (- 1/2 log(2 pi sigma_2^2) -(x-mu_2)^2/(2sigma_2^2)) p(x) dif x\
  &= - 1/2 log(2 pi sigma_2^2) integral_RR p(x) dif x - 1/(2sigma_2^2) integral_RR (x-mu_2)^2 p(x) dif x\
  &= - 1/2 log(2 pi sigma_2^2) - 1/(2sigma_2^2) [integral_RR x^2 p(x)-2mu_2 x p(x)+mu_2^2 p(x) dif x]\
$
We now use the notation $EE_p [g(x)] = integral g(x) p(x) dif x$, the expected
value of $g(x)$ under $p(x)$. We note that $EE_p [x^2] - EE_p [x]^2 =
"Var"(x)$, and thus $EE_p [x^2] = sigma_1^2 + mu_1^2$. Rewriting the last expression gives
$
  integral_RR p(x) log p(x) dif x & = - 1/2 log(2 pi sigma_2^2) - 1/(2sigma_2^2) [integral_RR x^2 p(x)-2mu_2 x p(x)-mu_2^2 p(x) dif x] \
  & = - 1/2 log(2 pi sigma_2^2) - 1/(2sigma_2^2) [EE_p [x^2] - 2mu_2 EE_p [x]+mu_2^2] \
  & = - 1/2 log(2 pi sigma_2^2) - 1/(2sigma_2^2) [sigma_1^2 + mu_1^2 - 2mu_2 mu_1 +mu_2^2] \
  & = - 1/2 log(2 pi sigma_2^2) - 1/(2sigma_2^2) [sigma_1^2 + (mu_1 - mu_2)^2] \
$
Putting everything together, we get
$
  dkl(P, Q) & = 1/2 log(sigma_2^2 / sigma_1^2) + sigma_1^2/(2sigma_2^2) + (mu_1 - mu_2)^2/(2sigma_2^2) - 1/2 \
            & =log(sigma_2 / sigma_1) + sigma_1^2/(2sigma_2^2) + (mu_1 - mu_2)^2/(2sigma_2^2) - 1/2. \
$
In particular, if $sigma_1 = sigma_2$, we have
$
  dkl(P, Q) & = (mu_1 - mu_2)^2/(2sigma_2^2) = dkl(Q, L).
$
We may also note that taking the limit as $sigma_1 -> 0^+$ (and fixing
$sigma_2$) yields an infinite value (since $1/x^2$ grows faster than $log x$ as
$x->0^+$), and the same can be said for $sigma_2 ->0^+$ (and fixing $sigma_1$).
This can be interpreted as it being equally as surprising if we obtain
information when none was emitted, as it is to not obtain any when it was sent.


#figure(
  image(width: 100%, "kl.png"),
  caption: [Level curves of the KL divergence for some values of $|mu_1 - mu_2|$. We note the *asymmetry* seen in the equation.],
)
#figure(
  image(width: 100%, "kl2.png"),
  caption: [Level curves of the KL divergence for some values of $sigma_1 = sigma_2 = sigma$. We note the *symmetry* seen in the equation.],
)

= Mutual Information
The *mutual information* of two random variables is a measure of the mutual
dependence between the two variables. It quantifies the "amount of information"
obtained about one random variable by observing the other random variable.

Let $(X,Y)$ be a pair of random variables with values over the space $cal(X)
times cal(Y)$. If their distribution is $P_(X,Y)$ and the marginal
distributions are $P_X$ and $P_Y$, the mutual information is defined as
$
  I(X;Y) = dkl(P_(X,Y), P_X times P_Y)
$
where $P_X times P_Y$ is the outer product distribution which assigns
probability $P_X (x) dot P_Y (y)$ to each $(x,y)$.

= Mutual Information of two Normal Distributions
Let $X = cal(N)(mu_X, sigma_X)$ and $Y = cal(N)(mu_Y, sigma_Y)$. Then, we recall the following:
$
  p(x,y) = 1/(2 pi sigma_X sigma_Y sqrt(1-rho^2)) exp(- 1/(1-rho^2)[(x-mu_X)^2/sigma_X^2 - 2 rho ((x - mu_X)(y-mu_Y))/(sigma_X sigma_Y) + (y - mu_Y)^2/sigma^2_Y])
$
We define $1/(2 pi sigma_X sigma_Y sqrt(1-rho^2)) =: kappa$ to shorten the expression.

Now, by definition of the KL divergence,
$
  I(X;Y) & = dkl(P_(X,Y), P_X times P_Y) \
  & = integral_(RR^2) p(x,y) log p(x,y)/(p(x)p(y)) dif x \
  & = integral_(RR^2) p(x,y) log p(x,y) dif x dif y - integral_(RR^2) p(x,y)log p(x) dif x dif y - integral_(RR^2) p(x,y) log p(y) dif x dif y \
  & = EE[log p(x,y)] - EE[log p(x)] - EE[log p(y)]
$
We then have that
$
  EE[log p(x,y)] &= log kappa - 1/(2 (1-rho^2)) [(EE[(x - mu_X)^2])/sigma_X^2 - (2rho)/(sigma_X sigma_Y) (rho sigma_X sigma_Y) + sigma_Y^2/sigma_X^2]\
  &=log kappa - 1/(2(1-rho^2))[1 - 2 rho^2 + 1]\
  &= log kappa - 1.
$
The other terms simplify in the following manner:
$
  integral_(RR) (integral_RR p(x,y) dif y) log p(x) dif x = integral_(RR) p(x) log p(x) dif x,
$
which we calculated earlier to be
$
  - 1/2 log(2 pi sigma_1^2) - 1/2.
$
Similarly,
$
  integral_(RR) (integral_RR p(x,y) dif x) log p(y) dif y = integral_(RR) p(y) log p(y) dif y
$
evaluates to
$
  - 1/2 log(2 pi sigma_2^2) - 1/2.
$
Altogether, we have
$
  I(X;Y) = log kappa + log (2 pi sigma_1 sigma_2) = log (2 pi sigma_X sigma_Y)/(2 pi sigma_X sigma_Y sqrt(1-rho^2)) = -1/2 log (1 - rho^2).
$
In particular, we notice that if $rho = 0$ (the variables have no correlation),
$
  I(X;Y) = 0.
$
If $rho -> 1^-$,
$
  I(X;Y) -> infinity.
$
And similarly, if $rho -> -1^+$,
$
  I(X;Y) -> infinity.
$
We also note that the mutual information does not depend on the means of each
distribution. This is due to the entropy only caring about the variance of the
data, and not the locality.
