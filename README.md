<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [['$', '$'] ],
 displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
 }
 });
</script>
# 1次元非定常移流拡散方程式の有限要素法によるソルバー
## 非定常移流拡散方程式とは
なんらかの物理量の一般的な流れを表す2階線形偏微分方程式であり，次の方程式で表される．

## オイラーの公式
オイラーの公式は以下のように与えられる。

$$ e^{i x} = \cos{x} + i \sin{x} $$

##  $ \varepsilon - \delta $ 論法
任意の $ \varepsilon > 0 $ についてある $ \delta > 0 $ が存在して、任意の $ x \in \mathbb{R} $ に対して $ 0 < |x - a| < \delta $ ならば $ |f(x) - f(a)| < \varepsilon $ を満たすとき $ f(x) $ は $ a $ で連続であるという。