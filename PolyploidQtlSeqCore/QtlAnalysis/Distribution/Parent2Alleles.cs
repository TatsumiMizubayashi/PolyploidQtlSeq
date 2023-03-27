namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 親2のアレル
    /// </summary>
    internal static class Parent2Alleles
    {
        /// <summary>
        /// 親2アレル配列を生成する。
        /// </summary>
        /// <param name="ploidy">倍数性</param>
        /// <param name="p2Plex">P2 Plex数</param>
        /// <returns>P2アレル配列</returns>
        public static Allele[] Generate(Ploidy ploidy, Parent2PlexNumber p2Plex)
        {
            if (ploidy.Value == p2Plex.Value)
                throw new ArgumentException($"{nameof(ploidy)} and {nameof(p2Plex)} have the same value.");

            var refCount = ploidy.Value - p2Plex.Value;
            var refAlleleQuery = Enumerable.Range(1, refCount)
                .Select(_ => new Allele(isAlt: false));
            var altAlleleQuery = Enumerable.Range(1, p2Plex.Value)
                .Select(_ => new Allele(isAlt: true));

            return refAlleleQuery.Concat(altAlleleQuery).ToArray();
        }
    }
}
