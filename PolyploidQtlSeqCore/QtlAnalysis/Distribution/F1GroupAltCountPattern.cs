namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// F1集団のAlt数パターン
    /// </summary>
    internal static class F1GroupAltCountPattern
    {
        /// <summary>
        /// F1集団のAlt数パターンを生成する。
        /// </summary>
        /// <param name="ploidy">倍数性</param>
        /// <param name="p2Plex">P2Plex数</param>
        /// <returns>F1集団Alt数パターン</returns>
        public static int[] Generate(Ploidy ploidy, Parent2PlexNumber p2Plex)
        {
            var p2Alleles = Parent2Alleles.Generate(ploidy, p2Plex);
            var k = ploidy.Value / 2;

            var f1AllelePatterns = Combination.Enumerate(p2Alleles, k).ToArray();

            var altNums = f1AllelePatterns
                .Select(alleles => alleles.Count(x => x.IsAlt))
                .OrderBy(x => x)
                .ToArray();

            return altNums;
        }
    }
}
