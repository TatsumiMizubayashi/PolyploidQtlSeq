namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 重複なしでN個の中からk個取り出す組み合わせ
    /// </summary>
    internal static class Combination
    {
        /// <summary>
        /// アレルからk個取り出す組み合わせを列挙する。
        /// </summary>
        /// <param name="alleles">アレル</param>
        /// <param name="k">取り出す数</param>
        /// <returns>組み合わせ</returns>
        public static IEnumerable<Allele[]> Enumerate(IEnumerable<Allele> alleles, int k)
        {
            if (k == 1)
            {
                foreach (var allele in alleles)
                {
                    yield return new[] { allele };
                }

                yield break;
            }

            foreach (var allele in alleles)
            {
                var leftside = new[] { allele };

                var unselectedAlleles = alleles
                    .SkipWhile(x => !ReferenceEquals(x, allele))
                    .Skip(1)    // allele自身をskip
                    .ToArray();

                foreach (var rightside in Enumerate(unselectedAlleles, k - 1))
                {
                    yield return leftside.Concat(rightside).ToArray();
                }
            }
        }
    }
}
