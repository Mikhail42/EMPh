package solution
object Solution {
  import math.exp
  import solution.Data

  val n0 = solution.Data.n0

  /** $e^{-k(k+1)\frac{K}{cR^2}t}$ */
  def expK(k: Int, t: Double) = exp(-(k*(k+1))*Data.kDelcR2*t);

  /** $P_{n}(x)=\frac{2n-1}{n}xP_{n-1}(x)-\frac{n-1}{n}P_{n-2}(x),\; \forall n>1$ */
  def getP_k(n: Int, pkm1: Double, pkm2: Double, x: Double) = ((n+n-1)*x*pkm1 - (n-1)*pkm2)/n;

  /** return result function $u(t,\theta)$
    * @param t --- time
    * @param theta --- angle, radians */
  def u(t: Double, theta: Double) = {
    val x = math.cos(theta)
    var pkm2 = 1.0;  // $P_{k-2}$
    var pkm1 = x;    // $P_{k-1}$
    var u    = 0.0;  // $u(t,\theta)$ --- result
    var pk   = 1.0;  // $P_{0}$
    for (k <- 2 to n0<<1) {
      pk =  getP_k(k,pkm1,pkm2,x)  // $P_{i}$
      pkm2 = pkm1
      pkm1 = pk
      // $\frac{(4k+1)P_{2k-2}(0)}{16k^2 (2k+1)(k+1)}\left(1 - e^{-2k(1+2k)\frac{K}{cR^2} t}\right)P_{2k}(\cos\theta)$
      if ((k & 1) == 0)
        u += ((k<<1)+1)*Data.P2k((k>>1)-1)*
          (1.0 - expK(k,t))*pk/
          (k.toLong*k*(k+1)*((k>>1)+1) << 2) //$2k->k$
    }
    u += 0.25*(1.0-expK(1,t))*x    //$\rfrac{1}{4}(1 - e^{-2\frac{K}{cR^2} t})\cos(\theta)$
    u *= Data.bR2DelK
    u += 0.25*Data.beta*t/Data.c   //$\frac{\beta}{4c} t$
    u
  }
}