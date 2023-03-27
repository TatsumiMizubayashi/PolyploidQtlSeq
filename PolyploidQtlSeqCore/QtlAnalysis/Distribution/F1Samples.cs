namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// F1集団のサンプル
    /// </summary>
    internal class F1Samples
    {
        /// <summary>
        /// F1集団サンプルを作成する。
        /// </summary>
        /// <param name="altCounts">各サンプルのAlt数</param>
        /// <param name="ploidy">ploidy</param>
        public F1Samples(int[] altCounts, Ploidy ploidy)
        {
            var altCount = altCounts.Sum();

            AltFrequency = (double)altCount / (altCounts.Length * ploidy.Value);
            RefFrequency = 1.0 - AltFrequency;
        }

        /// <summary>
        /// Alt割合を取得する。
        /// </summary>
        public double AltFrequency { get; }

        /// <summary>
        /// Ref割合を取得する。
        /// </summary>
        public double RefFrequency { get; }

        /// <summary>
        /// リード配列に変換する。
        /// </summary>
        /// <param name="depth">Depth</param>
        /// <returns>Reads</returns>
        public Reads ToReads(int depth)
        {
            var alleles = new[] { 0, 1 };
            var alleleProbability = new[] { RefFrequency, AltFrequency };

            var sampling = new RandomSampling();
            var gtValues = sampling.Choice(alleles, alleleProbability, depth);

            return new Reads(gtValues);
        }
    }
}
