package solution
object Data{
  var c: Double = 2.0
  var K: Double = 0.01
  var R: Double = 2.0
  var beta: Double = 0.1
  // @param eps: accuracy: $\abs{u-u^*}, u^* = \sum_{0}^{N} \{\dots\}$
  var eps: Double = 1e-7

  def R2 = R*R
  def kDelcR2 = K/(c*R2)
  def bR2DelK = beta*R2/K

  // Long.Max = 2^63. I can use n0^4,  => n0 < (2^63)^0.25 = 55108.
  val maxN0 = 40000
  /** counts of $\{f_{2k}\}$, $\{P_{2k}(0)\}$ */
  /** $n_0 = \Bigg\lceil \sqrt{\frac{\beta R^2}{16K\varepsilon}}\Bigg\rceil$*/
  def n0 =  math.min(math.sqrt(bR2DelK/(16*eps)), (maxN0>>1)).toInt

  /** $P_{2n}(0) = - P_{2n-2}(0) \frac{2n-1}{2n},\; P_{0}(0)=1$ */
  val P2k = new Array[Double](n0)
  def setAllP2k{
    P2k(0)=1
    for (i <- 1 until n0)
      P2k(i) = -P2k(i-1)*((i << 1)-1)/(i << 1)
  }

  // time
  var t: Double = 1e-3
  // initial of $\{P_{2k}\}$
  def init = setAllP2k
}
