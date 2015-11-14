package solution
object Solution {
  import math.exp
  import solution.Data

  /** $e^{-k(k+1)\frac{K}{cR^2}t}$ */
  def expK(k: Int, t: Double, kDelcR2: Double) = exp(-(k*(k+1))*kDelcR2*t);

  /** $P_{n}(x)=\frac{2n-1}{n}xP_{n-1}(x)-\frac{n-1}{n}P_{n-2}(x),\; \forall n>1$ */
  def getP_k(n: Int, pkm1: Double, pkm2: Double, x: Double) = ((n+n-1)*x*pkm1 - (n-1)*pkm2)/n;

  /**
    * computes of u(t,x) on [0,2pi)
    * @param t --- time
    * @param countPoints --- the number of points on [0,2pi]
    * @return (0:2pi/n:2pi, u(t,0:2pi/n:2pi)), where n = countPoint
    */
  def masU(t: Double, countPoints: Int) = {
    val n       = countPoints
    val n0      = Data.n0
    val kDelcR2 = Data.kDelcR2
    val bR2DelK = Data.bR2DelK
    val beta    = Data.beta
    val c       = Data.c

    /** return result function $u(t,\theta)$
      * @param theta --- angle, radians */
    def u(theta: Double) = {
      val x       = math.cos(theta)
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
            (1.0 - expK(k,t,kDelcR2))*pk/
            (k.toLong*k*(k+1)*((k>>1)+1) << 2) //$2k->k$
      }

      u += 0.25*(1.0-expK(1,t,kDelcR2))*x    //$\rfrac{1}{4}(1 - e^{-2\frac{K}{cR^2} t})\cos(\theta)$
      u *= bR2DelK
      u += 0.25*beta*t/c            //$\frac{\beta}{4c} t$
      u
    }

    val X = new Array[Double](n)
    val Y = new Array[Double](n)
    X(0) = 0.0; Y(0) = u(X(0))
    val h: Double = math.Pi/n
    for (i <- 1 until n) {
      X(i) = X(i-1)+h
      Y(i) = u(X(i))
    }

    (X,Y)
  }
}
